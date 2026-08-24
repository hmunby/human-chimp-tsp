"""
BMAP-matched null-locus background generator for the locus -> gene enrichment.

Rationale: the test loci are ascertained in high-BMAP (low background-selection)
regions, and the gene-assignment step (L2G -> supported-molQTL -> nearest) has its
own structural biases. To make a valid DAVID background we resample null loci that
MATCH the foreground on BMAP + locus span + SNP count, then push them through the
IDENTICAL assignment rule (assign_one_gene) and collect the genes.

This script has two stages:
  STAGE 1 (this file, no network): sample null loci, validate BMAP matching,
          write a manifest. Run and review BEFORE the expensive API stage.
  STAGE 2 (added after sign-off): pull ref/alt from the VCF for sampled SNPs,
          query OpenTargets (cached), assign one gene per null locus, emit the
          background gene universe + per-replicate null gene sets.

Foreground loci: ~500 kb windows, 1-17 SNPs, median BMAP 936.
"""

import argparse
import json
import subprocess
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd

# Reuse the locked assignment rule + OT plumbing from the foreground pipeline
import opentargets_assign as ota

# ── Paths ─────────────────────────────────────────────────────────────────────

HERE = Path(__file__).resolve().parent
GA_DIR = HERE                       # null caches and gene-set outputs are written here

# The reported foreground, from 4_TSPRegions/3_GeneAssignment/.
_FG = HERE / "../../../4_SINGER_GenealogicalReconstruction_TSP/4_TSPRegions/3_GeneAssignment"
REGION_TSV  = _FG / "region_gene_assignments.pooled.tsv"      # foreground loci
VARIANT_TSV = _FG / "variant_gene_assignments.pooled.tsv"     # foreground per-SNP BMAP

# All common autosomal SNPs with their BMAP (chrom, start, end, BMAP); not bundled.
POOL_BED = HERE / "resources/snps_autosomes.bed"
MANIFEST_TSV = GA_DIR / "null_loci_manifest.tsv"
NULL_CACHE   = GA_DIR / "ot_null_cache.json"        # OT responses for null SNPs
NULL_ASSIGN_TSV = GA_DIR / "null_locus_assignments.tsv"
NULL_SETS_TSV   = GA_DIR / "null_gene_sets.tsv"     # one gene-set per replicate
CHECKPOINT_EVERY = 100   # batches between full cache writes (see fetch_null_ot)
BACKGROUND_TXT  = GA_DIR / "background_genes.txt"   # DAVID custom background universe

VCF = "resources/1000GP.wgs.biallelic.sites.YRI.CpG.snpeff.MR.hg38.MAF05.vcf.gz"

CHR_LENGTHS = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468,
}

BMAP_TOL = 20   # anchor BMAP must be within +/- this of the foreground locus mean


# ── Load foreground ───────────────────────────────────────────────────────────

def load_foreground():
    fg = pd.read_csv(REGION_TSV, sep="\t")
    # ACTUAL SNP-cluster span (not the 500 kb padding): max-min of the SNP positions
    v = pd.read_csv(VARIANT_TSV, sep="\t")[["region", "POSITION_HG38", "BMAP"]]
    v["POSITION_HG38"] = v["POSITION_HG38"].astype(int)
    v["BMAP"] = pd.to_numeric(v["BMAP"], errors="coerce")
    span = (v.groupby("region")["POSITION_HG38"]
              .agg(lambda s: int(s.max() - s.min())).rename("span"))
    fg = fg.merge(span, on="region", how="left")
    fg["span"] = fg["span"].fillna(0).astype(int)
    return fg, v


# ── Load genome-wide common-SNP pool, indexed per chromosome ──────────────────

def load_pool():
    print(f"Loading SNP pool: {POOL_BED}")
    pool = pd.read_csv(POOL_BED, sep="\t", header=None,
                       names=["chrom", "start", "pos", "bmap"],
                       dtype={"chrom": str, "pos": np.int64, "bmap": np.float32},
                       usecols=["chrom", "pos", "bmap"])
    pool = pool[pool["chrom"].isin(CHR_LENGTHS)].reset_index(drop=True)
    # per-chrom sorted arrays for fast window queries
    by_chrom = {}
    for c, g in pool.groupby("chrom"):
        gp = g.sort_values("pos")
        by_chrom[c] = (gp["pos"].to_numpy(), gp["bmap"].to_numpy())
    print(f"  pool: {len(pool):,} SNPs across {len(by_chrom)} chromosomes")
    return pool, by_chrom


# ── Sampling ──────────────────────────────────────────────────────────────────

def sample_null_loci(fg, by_chrom, pool, n_rep, seed=0):
    """For each replicate, build one matched null locus per foreground locus.

    A null locus = a BMAP-matched anchor SNP + its (n-1) nearest common-SNP
    neighbours (a tight cluster). The actual SNP-cluster span is the target window;
    if too few common SNPs fall inside it, the window expands to the nearest n.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for _, locus in fg.iterrows():
        target_bmap = locus["mean_bmap"]
        span        = int(locus["span"])
        n_snps      = int(locus["n_snps"])
        half        = max(span // 2, 1)

        # BMAP-matched candidate anchors
        cand = pool[np.abs(pool["bmap"] - target_bmap) <= BMAP_TOL]
        if len(cand) == 0:
            cand = pool.iloc[(pool["bmap"] - target_bmap).abs().argsort()[:100000]]
        anchors = cand.sample(n=n_rep, replace=True,
                              random_state=int(rng.integers(1_000_000_000)))

        for rep, (_, anc) in enumerate(anchors.iterrows()):
            c = anc["chrom"]
            apos = int(anc["pos"])
            positions, bmaps = by_chrom[c]
            ai = int(np.searchsorted(positions, apos))

            if n_snps == 1:
                sel = np.array([min(ai, len(positions) - 1)])
            else:
                lo = np.searchsorted(positions, apos - half, "left")
                hi = np.searchsorted(positions, apos + half, "right")
                if hi - lo >= n_snps:                 # enough within the cluster span
                    win = np.arange(lo, hi)
                else:                                  # expand to nearest neighbours
                    lo = max(0, ai - n_snps * 10)
                    hi = min(len(positions), ai + n_snps * 10)
                    win = np.arange(lo, hi)
                d = np.abs(positions[win] - apos)
                sel = win[np.argsort(d)[:n_snps]]      # the n nearest to the anchor

            sel_pos  = np.sort(positions[sel])
            sel_bmap = bmaps[sel]
            rows.append(dict(
                replicate=rep,
                fg_region=locus["region"],
                fg_gene=locus["assigned_gene"],
                fg_tier=locus["assignment_tier"],
                target_bmap=target_bmap, n_snps=n_snps, target_span=span,
                null_chrom=c, anchor_pos=apos,
                realized_span=int(sel_pos.max() - sel_pos.min()),
                null_positions=",".join(map(str, sel_pos.tolist())),
                null_bmap_mean=round(float(sel_bmap.mean()), 1),
                null_bmaps=",".join(str(int(b)) for b in sel_bmap),
            ))
    return pd.DataFrame(rows)


# ── Validation ────────────────────────────────────────────────────────────────

def validate(manifest, fg_snp_bmap):
    fg_b = fg_snp_bmap["BMAP"].dropna().to_numpy()
    null_b = np.concatenate([
        np.array(s.split(","), dtype=float) for s in manifest["null_bmaps"]
    ])
    print("\n── BMAP match (per-SNP) ──")
    qs = [0, 10, 25, 50, 75, 90, 100]
    print("  pct :   " + "  ".join(f"{q:>5}" for q in qs))
    print("  fore:   " + "  ".join(f"{np.percentile(fg_b,q):5.0f}" for q in qs))
    print("  null:   " + "  ".join(f"{np.percentile(null_b,q):5.0f}" for q in qs))
    print(f"  foreground SNPs: {len(fg_b)}  |  null SNPs: {len(null_b):,}")
    print(f"  mean  fore {fg_b.mean():.1f}  vs  null {null_b.mean():.1f}")
    from scipy import stats
    ks = stats.ks_2samp(fg_b, null_b)
    print(f"  KS test fore vs null: D={ks.statistic:.3f}, p={ks.pvalue:.3f} "
          f"(high p = well matched)")

    print("\n── Cluster span match (multi-SNP loci) ──")
    multi = manifest[manifest["n_snps"] > 1]
    print(f"  target span  : median {multi['target_span'].median():.0f} bp")
    print(f"  realized span: median {multi['realized_span'].median():.0f} bp")
    over = (multi["realized_span"] > multi["target_span"] * 2).mean()
    print(f"  null loci forced >2x wider than target (sparse common SNPs): {over*100:.0f}%")


# ── Stage 2: alleles -> OpenTargets -> assign one gene per null locus ──────────

def fetch_alleles(uniq):
    """uniq: DataFrame[chrom,pos]; return {(chrom,pos): (ref,alt)} from the VCF."""
    print(f"  looking up ref/alt for {len(uniq):,} unique null SNPs from VCF…")
    with tempfile.TemporaryDirectory() as tmp:
        bed = Path(tmp) / "pos.bed"
        out = Path(tmp) / "alleles.txt"
        b = uniq.copy()
        b["start"] = b["pos"].astype(int) - 1
        b[["chrom", "start", "pos"]].sort_values(["chrom", "start"]).to_csv(
            bed, sep="\t", header=False, index=False)
        subprocess.run(
            f"bcftools query -R {bed} -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' {VCF} > {out}",
            shell=True, check=True)
        al = pd.read_csv(out, sep="\t", header=None,
                         names=["chrom", "pos", "ref", "alt"], dtype=str)
    al["alt"] = al["alt"].str.split(",").str[0]          # biallelic: single ALT
    return {(r.chrom, int(r.pos)): (r.ref, r.alt) for r in al.itertuples()}


def fetch_null_ot(variant_ids):
    """Query OT for null variant IDs with the same batching, into a SEPARATE cache."""
    cache = json.loads(NULL_CACHE.read_text()) if NULL_CACHE.exists() else {}
    todo = [v for v in variant_ids if v not in cache]
    print(f"  OT: {len(todo):,} to fetch ({len(variant_ids)-len(todo):,} cached)")
    for i in range(0, len(todo), ota.BATCH_SIZE):
        chunk = todo[i:i + ota.BATCH_SIZE]
        aliases = {f"v{j}": v for j, v in enumerate(chunk)}
        q = "{" + "".join(ota._variant_subquery(a, v) for a, v in aliases.items()) + "}"
        data = ota._post(q)
        for a, v in aliases.items():
            cache[v] = data.get(a)
        # Checkpoint PERIODICALLY, not after every batch. Rewriting the whole cache each batch is
        # O(n) per batch and therefore quadratic overall: at 2,000 replicates the cache reaches
        # ~574 MB and per-batch writes would mean ~3.2 TB of I/O and ~32 h, against ~0.9 h of
        # actual network time. Every 100 batches costs ~32 GB total and risks losing at most
        # 100*BATCH_SIZE variants on a crash, which are simply re-fetched on resume.
        nb = i // ota.BATCH_SIZE
        if nb % CHECKPOINT_EVERY == 0:
            NULL_CACHE.write_text(json.dumps(cache))
        if nb % 100 == 0:
            print(f"    {min(i+ota.BATCH_SIZE, len(todo)):,}/{len(todo):,}", flush=True)
        time.sleep(0.3)
    NULL_CACHE.write_text(json.dumps(cache))   # final flush
    return cache


def stage2_assign(manifest):
    """Assign one gene to each null locus using the locked foreground rule."""
    # explode null loci into (locus_uid, chrom, pos)
    recs = []
    for row in manifest.itertuples():
        uid = f"{row.fg_region}|rep{row.replicate}"
        for p in str(row.null_positions).split(","):
            recs.append((uid, row.fg_region, row.replicate, row.null_chrom, int(p)))
    snps = pd.DataFrame(recs, columns=["locus_uid", "fg_region", "replicate", "chrom", "pos"])

    uniq = snps[["chrom", "pos"]].drop_duplicates()
    alleles = fetch_alleles(uniq)

    def vid(r):
        a = alleles.get((r.chrom, r.pos))
        # OpenTargets variant IDs have NO 'chr' prefix (e.g. 15_66628096_A_G)
        c = r.chrom.replace("chr", "")
        return f"{c}_{r.pos}_{a[0]}_{a[1]}" if a else None
    snps["variant_id"] = [vid(r) for r in snps.itertuples()]
    snps = snps.dropna(subset=["variant_id"])

    # OT functional evidence (cached) for all unique null variants
    cache = fetch_null_ot(snps["variant_id"].drop_duplicates().tolist())

    # nearest gene (one bedtools run; reuse foreground's GENCODE-PC procedure)
    print("  nearest gene (bedtools closest vs GENCODE v49)…")
    snp_df = snps.drop_duplicates("variant_id").copy()
    snp_df["CHROM_HG38"]   = snp_df["chrom"].str.replace("chr", "", regex=False)
    snp_df["POSITION_HG38"] = snp_df["pos"]
    near = ota.nearest_genes(snp_df[["CHROM_HG38", "POSITION_HG38", "variant_id"]])
    near_idx = near.set_index("variant_id")

    # one gene per null locus
    print(f"  assigning {snps['locus_uid'].nunique():,} null loci…")
    rows = []
    for uid, g in snps.groupby("locus_uid"):
        vids = g["variant_id"].tolist()
        near_sub = near_idx.reindex(vids).reset_index()
        a = ota.assign_one_gene(vids, cache, near_sub)
        rows.append(dict(locus_uid=uid, fg_region=g["fg_region"].iloc[0],
                         replicate=int(g["replicate"].iloc[0]),
                         null_chrom=g["chrom"].iloc[0], n_snps=len(g),
                         assigned_gene=a["assigned_gene"], assignment_tier=a["assignment_tier"]))
    null_tbl = pd.DataFrame(rows)
    null_tbl.to_csv(NULL_ASSIGN_TSV, sep="\t", index=False)

    # per-replicate gene sets (for empirical-null option) + background universe
    sets = (null_tbl.groupby("replicate")["assigned_gene"]
                    .apply(lambda s: ",".join(sorted(set(s.dropna())))).reset_index())
    sets.to_csv(NULL_SETS_TSV, sep="\t", index=False)

    # DAVID requires foreground ⊆ background -> union the foreground genes in
    fg_genes = set(pd.read_csv(REGION_TSV, sep="\t")["assigned_gene"].dropna())
    null_genes = set(null_tbl["assigned_gene"].dropna())
    background = sorted(null_genes | fg_genes)
    BACKGROUND_TXT.write_text("\n".join(background) + "\n")
    print(f"  background = {len(null_genes)} null genes + "
          f"{len(fg_genes - null_genes)} foreground-only genes = {len(background)}")

    # ── report ────────────────────────────────────────────────────────────────
    print(f"\nSaved: {NULL_ASSIGN_TSV}")
    print(f"Saved: {NULL_SETS_TSV}")
    print(f"Saved DAVID background ({len(background)} genes): {BACKGROUND_TXT}")
    print(f"\n  null-locus tiers: {null_tbl['assignment_tier'].value_counts().to_dict()}")
    fg = pd.read_csv(REGION_TSV, sep="\t")
    fg_tiers = fg["assignment_tier"].value_counts(normalize=True).mul(100).round(0).to_dict()
    nl_tiers = null_tbl["assignment_tier"].value_counts(normalize=True).mul(100).round(0).to_dict()
    print(f"  tier % foreground: {fg_tiers}")
    print(f"  tier % null      : {nl_tiers}")
    return null_tbl, background


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-rep", type=int, default=100,
                    help="null replicates per foreground locus")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--assign", action="store_true",
                    help="Stage 2: hit OpenTargets + assign genes (slow, external)")
    args = ap.parse_args()

    fg, fg_snp_bmap = load_foreground()
    pool, by_chrom = load_pool()

    print(f"\nSampling {args.n_rep} null replicates x {len(fg)} loci "
          f"= {args.n_rep * len(fg):,} null loci…")
    manifest = sample_null_loci(fg, by_chrom, pool, args.n_rep, args.seed)

    n_null_snps = manifest["null_positions"].str.count(",").add(1).sum()
    print(f"  built {len(manifest):,} null loci  ({n_null_snps:,} null SNPs total)")

    manifest.to_csv(MANIFEST_TSV, sep="\t", index=False)
    print(f"  saved manifest: {MANIFEST_TSV}")

    validate(manifest, fg_snp_bmap)

    if args.assign:
        print("\n=== STAGE 2: allele lookup -> OpenTargets -> gene assignment ===")
        stage2_assign(manifest)
    else:
        print("\n(Stage 1 only. Re-run with --assign to query OpenTargets and "
              "build the background.)")


if __name__ == "__main__":
    main()
