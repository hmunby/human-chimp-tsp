"""
Three-tier variant -> gene assignment, reproducing the manual OpenTargets workflow.

Tier 1 `l2g`      GWAS locus-to-gene effector prediction scoring >= L2G_THRESH. L2G already
                  integrates molQTL colocalisation, distance, chromatin and VEP, so it takes
                  precedence over the individual lines of evidence below.
Tier 2 `molqtl`   membership of a molecular-QTL (eQTL/sQTL/pQTL/...) 95% credible set, assigned
                  to the QTL-regulated gene. Where a region has several molQTL genes the winner
                  is chosen on breadth: tissues, then modalities, then credible sets.
Tier 3 `nearest`  bedtools closest against GENCODE v49 protein-coding genes (autosomes) -- the
                  same proximity rule as the glycoprotein/territory analysis.

The rule is applied per REGION on evidence aggregated across its TSPs, and each variant inherits
its region's decision, so the per-variant and per-region tables cannot disagree.

Each variant's posterior probability (PIP) within every credible set and the credible set's
confidence are recorded, so a PIP threshold can be applied later.

Design notes verified against the live API (2026-06):
  * OpenTargets Genetics was merged into the Platform; endpoint is
    https://api.platform.opentargets.org/api/v4/graphql
  * variant IDs are CHR_POS_REF_ALT on GRCh38 (e.g. 1_154453788_C_T)
  * CredibleSet.studyType / .qtlGeneId give QTL type + regulated gene directly
  * CredibleSet.locus(variantIds:[id]) -> Locus.posteriorProbability = this
    variant's PIP within that set
  * ~40% of arbitrary SNPs are not in the OT variant index (return null) -> these
    legitimately have no functional evidence and fall through to the nearest-gene tier.

Outputs (`.pooled` tag shown; --snp-table and --out-* select an alternative set):
  ot_variant_cache.json                 - raw API responses keyed by variant_id (reusable)
  variant_gene_assignments.pooled.tsv   - tidy per-variant table (117 TSPs)
  region_gene_assignments.pooled.tsv    - one gene per TSP region (59), the foreground set
"""

import json
import time
import urllib.request
import urllib.error
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

# ---- inputs (edit) -------------------------------------------------------------------------
HERE = Path(__file__).resolve().parent

# The TSPs with their cluster region ids, from ../2_ClusteringLD/aggregate_seeds.py
SNP_TABLE = HERE / "../2_ClusteringLD/tsps.pooled.for_gene_assignment.txt"

# GENCODE v49 protein-coding genes, autosomes, BED with gene names (Tier 3 nearest gene).
# Not bundled -- see README.
GENCODE_BED = HERE / "resources/gencode.v49.annotation.genes.protein_coding.autosomes.gene_names.bed"

OUT_DIR = HERE
CACHE_JSON = OUT_DIR / "ot_variant_cache.json"          # raw API responses, reusable
OUT_TSV = OUT_DIR / "variant_gene_assignments.pooled.tsv"      # per-variant (transparency)
REGION_TSV = OUT_DIR / "region_gene_assignments.pooled.tsv"    # per-locus, ONE gene (foreground)

OT_API = "https://api.platform.opentargets.org/api/v4/graphql"

# molecular-QTL study types treated as functional evidence (tier 2)
MOLQTL_TYPES = {"eqtl", "pqtl", "sqtl", "tuqtl", "edqtl", "sceqtl", "scsqtl"}

BATCH_SIZE = 20   # variants per aliased GraphQL POST

L2G_THRESH = 0.5    # min L2G score to accept a GWAS effector-gene prediction (tier 1)


# ── GraphQL ────────────────────────────────────────────────────────────────────

def _variant_subquery(alias: str, vid: str) -> str:
    """One aliased variant block; locus filtered to THIS variant for its PIP."""
    return f'''
    {alias}: variant(variantId: "{vid}") {{
      id
      rsIds
      mostSevereConsequence {{ label }}
      credibleSets(page: {{index: 0, size: 200}}) {{
        count
        rows {{
          studyLocusId
          studyType
          qtlGeneId
          confidence
          study {{ id studyType target {{ id approvedSymbol }} biosample {{ biosampleId biosampleName }} }}
          locus(variantIds: ["{vid}"], page: {{index: 0, size: 1}}) {{
            rows {{ posteriorProbability is95CredibleSet }}
          }}
          l2GPredictions(page: {{index: 0, size: 3}}) {{
            rows {{ score target {{ id approvedSymbol }} }}
          }}
        }}
      }}
    }}'''


def _post(query: str, retries: int = 4) -> dict:
    body = json.dumps({"query": query}).encode()
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                OT_API, data=body,
                headers={"Content-Type": "application/json"},
            )
            with urllib.request.urlopen(req, timeout=90) as r:
                payload = json.load(r)
            if "errors" in payload:
                # validation errors won't fix on retry; surface immediately
                raise RuntimeError(f"GraphQL errors: {payload['errors'][:2]}")
            return payload["data"]
        except (urllib.error.URLError, TimeoutError) as e:
            last = e
            time.sleep(2 ** attempt)  # backoff: 1,2,4,8s
    raise RuntimeError(f"OpenTargets request failed after {retries} tries: {last}")


def fetch_variants(variant_ids):
    """Query OT for all variant_ids with caching + alias batching."""
    cache = {}
    if CACHE_JSON.exists():
        cache = json.loads(CACHE_JSON.read_text())
        print(f"  Loaded {len(cache):,} cached variants from {CACHE_JSON.name}")

    todo = [v for v in variant_ids if v not in cache]
    print(f"  {len(todo):,} variants to fetch ({len(variant_ids) - len(todo):,} cached)")

    for i in range(0, len(todo), BATCH_SIZE):
        chunk = todo[i:i + BATCH_SIZE]
        aliases = {f"v{j}": vid for j, vid in enumerate(chunk)}
        query = "{" + "".join(_variant_subquery(a, v) for a, v in aliases.items()) + "}"
        data = _post(query)
        for a, vid in aliases.items():
            cache[vid] = data.get(a)   # None if not in OT index
        CACHE_JSON.write_text(json.dumps(cache))
        print(f"  fetched {min(i + BATCH_SIZE, len(todo)):,}/{len(todo):,}")
        time.sleep(0.3)               # be polite to the public API

    return cache


# ── Parse OT response into functional-evidence rows ───────────────────────────

def parse_variant(vid: str, rec: dict) -> dict:
    """Collapse one variant's OT record into assignment-relevant fields."""
    out = dict(
        variant_id=vid, rsid="", in_ot_index=(rec is not None),
        most_severe="", n_credible_sets=0,
        molqtl_genes="", molqtl_details="",      # tier-2 evidence
        gwas_l2g_genes="",                        # tier-1 evidence
        best_molqtl_gene="", best_molqtl_pip="", best_molqtl_type="", best_molqtl_tissue="",
    )
    if rec is None:
        return out

    out["rsid"] = ",".join(rec.get("rsIds") or [])
    msc = rec.get("mostSevereConsequence") or {}
    out["most_severe"] = msc.get("label", "")
    cs = rec.get("credibleSets") or {}
    rows = cs.get("rows") or []
    out["n_credible_sets"] = cs.get("count", len(rows))

    molqtl = []     # (gene_symbol, gene_id, qtl_type, tissue, pip, confidence)
    gwas_l2g = []   # (gene_symbol, l2g_score)

    for r in rows:
        stype = (r.get("studyType") or "").lower()
        study = r.get("study") or {}
        tgt   = study.get("target") or {}
        bios  = study.get("biosample") or {}
        loc_rows = ((r.get("locus") or {}).get("rows") or [])
        pip = loc_rows[0].get("posteriorProbability") if loc_rows else None

        if stype in MOLQTL_TYPES:
            sym = tgt.get("approvedSymbol") or r.get("qtlGeneId") or ""
            molqtl.append((
                sym, tgt.get("id") or r.get("qtlGeneId") or "",
                stype, bios.get("biosampleName") or "",
                pip, r.get("confidence") or "",
            ))
        elif stype == "gwas":
            for p in ((r.get("l2GPredictions") or {}).get("rows") or []):
                pt = p.get("target") or {}
                if pt.get("approvedSymbol"):
                    gwas_l2g.append((pt["approvedSymbol"], p.get("score")))

    # de-dup molQTL genes, keep highest-PIP credible set per gene.
    # stored tuple: (sym, gid, stype, tissue, pip_or_None, conf, sort_score)
    by_gene = {}
    for sym, gid, stype, tissue, pip, conf in molqtl:
        key = sym or gid
        if not key:
            continue
        score = pip if pip is not None else -1.0
        prev = by_gene.get(key)
        if prev is None or score > prev[6]:
            by_gene[key] = (sym, gid, stype, tissue, pip, conf, score)

    if by_gene:
        out["molqtl_genes"] = ",".join(sorted(by_gene))
        out["molqtl_details"] = ";".join(
            f"{s}|{t}|{tis}|pip={'NA' if pp is None else format(pp, '.3g')}|{c}"
            for (s, gid, t, tis, pp, c, sc) in by_gene.values()
        )
        # best = highest PIP across genes
        best = max(by_gene.values(), key=lambda x: x[6])
        out["best_molqtl_gene"]   = best[0]
        out["best_molqtl_pip"]    = "" if best[4] is None else round(best[4], 4)
        out["best_molqtl_type"]   = best[2]
        out["best_molqtl_tissue"] = best[3]

    # de-dup GWAS L2G genes, keep best score
    g_by_gene = {}
    for sym, sc in gwas_l2g:
        sc = sc if sc is not None else -1
        if sym not in g_by_gene or sc > g_by_gene[sym]:
            g_by_gene[sym] = sc
    if g_by_gene:
        out["gwas_l2g_genes"] = ",".join(
            f"{s}({sc:.2f})" for s, sc in sorted(g_by_gene.items(), key=lambda x: -x[1])
        )

    return out


# ── Tier 3: bedtools closest vs GENCODE v49 ───────────────────────────────────

CHR_LENGTHS = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468,
}


def nearest_genes(df: pd.DataFrame) -> pd.DataFrame:
    """Return variant_id -> (nearest_gene, nearest_distance) via bedtools closest."""
    with tempfile.TemporaryDirectory() as tmp:
        snp_bed   = Path(tmp) / "snps.bed"
        gen_srt   = Path(tmp) / "gencode_sorted.bed"
        out_bed   = Path(tmp) / "closest.bed"
        genome_f  = Path(tmp) / "genome.txt"

        # genome file fixes bedtools closest when the SNP set spans fewer
        # chromosomes than GENCODE (defines a consistent global order)
        with open(genome_f, "w") as fh:
            for i in range(1, 23):
                fh.write(f"chr{i}\t{CHR_LENGTHS[f'chr{i}']}\n")

        bed = df[["CHROM_HG38", "POSITION_HG38", "variant_id"]].copy()
        bed["chrom"] = "chr" + bed["CHROM_HG38"].astype(str)
        bed["start"] = bed["POSITION_HG38"].astype(int) - 1
        bed["end"]   = bed["POSITION_HG38"].astype(int)
        bed[["chrom", "start", "end", "variant_id"]].sort_values(
            ["chrom", "start"]
        ).to_csv(snp_bed, sep="\t", header=False, index=False)

        subprocess.run(f"sort -k1,1V -k2,2n {GENCODE_BED} > {gen_srt}",
                       shell=True, check=True)
        subprocess.run(
            f"sort -k1,1V -k2,2n {snp_bed} | "
            f"bedtools closest -a stdin -b {gen_srt} -d -g {genome_f} > {out_bed}",
            shell=True, check=True,
        )
        closest = pd.read_csv(
            out_bed, sep="\t", header=None,
            names=["chrom", "start", "end", "variant_id",
                   "g_chrom", "g_start", "g_end", "nearest_gene", "nearest_distance"],
        )
    # a SNP can tie two genes -> keep the first
    closest = closest.drop_duplicates("variant_id")
    return closest[["variant_id", "nearest_gene", "nearest_distance"]]


# ── Region-level evidence + one-gene-per-locus assignment ─────────────────────

def region_evidence(variant_ids, cache):
    """Aggregate molQTL multiplicity and GWAS-L2G across a locus's variants.

    Returns:
      molq: {gene: {'sets':int, 'tissues':set, 'mods':set}}
      l2g : {gene: max L2G score across the region's GWAS credible sets}
    """
    from collections import defaultdict
    molq = defaultdict(lambda: {"sets": 0, "tissues": set(), "mods": set()})
    l2g  = {}
    for vid in variant_ids:
        rec = cache.get(vid)
        if not rec:
            continue
        for r in (rec.get("credibleSets") or {}).get("rows") or []:
            st    = (r.get("studyType") or "").lower()
            study = r.get("study") or {}
            if st in MOLQTL_TYPES:
                tgt  = study.get("target") or {}
                bios = study.get("biosample") or {}
                sym  = tgt.get("approvedSymbol") or r.get("qtlGeneId")
                if not sym:
                    continue
                molq[sym]["sets"]    += 1
                molq[sym]["tissues"].add(bios.get("biosampleName") or "?")
                molq[sym]["mods"].add(st)
            elif st == "gwas":
                for p in ((r.get("l2GPredictions") or {}).get("rows") or []):
                    pt = p.get("target") or {}
                    sym, sc = pt.get("approvedSymbol"), p.get("score")
                    if sym and sc is not None:
                        l2g[sym] = max(l2g.get(sym, 0.0), sc)
    return molq, l2g


def assign_one_gene(variant_ids, cache, near_sub):
    """Apply L2G -> supported-molQTL -> nearest, returning a single gene per locus."""
    molq, l2g = region_evidence(variant_ids, cache)

    # nearest gene for the locus: majority across SNPs, tie-broken by min distance
    nn = near_sub.dropna(subset=["nearest_gene"])
    if len(nn):
        agg = (nn.groupby("nearest_gene")
                 .agg(n=("nearest_gene", "size"), mind=("nearest_distance", "min"))
                 .reset_index()
                 .sort_values(["n", "mind"], ascending=[False, True]))
        nearest_gene = agg.iloc[0]["nearest_gene"]
        nearest_dist = int(nn[nn["nearest_gene"] == nearest_gene]["nearest_distance"].min())
    else:
        nearest_gene, nearest_dist = "", -1

    # multiplicity string for ALL molQTL genes (recorded regardless of tier)
    mult_str = ";".join(
        f"{s}|sets={i['sets']}|tissues={len(i['tissues'])}|mods={'+'.join(sorted(i['mods']))}"
        for s, i in sorted(molq.items(), key=lambda x: -len(x[1]["tissues"]))
    )
    l2g_str = ",".join(f"{s}({v:.2f})" for s, v in sorted(l2g.items(), key=lambda x: -x[1]))

    # Tier 1: confident GWAS L2G
    l2g_hi = {k: v for k, v in l2g.items() if v >= L2G_THRESH}
    if l2g_hi:
        gene = max(l2g_hi, key=l2g_hi.get)
        return dict(assigned_gene=gene, assignment_tier="l2g",
                    assignment_evidence=f"L2G={l2g_hi[gene]:.2f}",
                    nearest_gene=nearest_gene, nearest_distance=nearest_dist,
                    n_molqtl_genes=len(molq), molqtl_multiplicity=mult_str, l2g_genes=l2g_str)

    # Tier 2: ANY molQTL gene beats the nearest gene. Multiplicity (tissues, then
    # modalities, then sets) is used ONLY to choose among multiple molQTL genes;
    # ties are broken alphabetically for determinism.
    if molq:
        gene = sorted(molq, key=lambda s: (-len(molq[s]["tissues"]),
                                           -len(molq[s]["mods"]),
                                           -molq[s]["sets"], s))[0]
        i = molq[gene]
        return dict(assigned_gene=gene, assignment_tier="molqtl",
                    assignment_evidence=f"sets={i['sets']},tissues={len(i['tissues'])},mods={'+'.join(sorted(i['mods']))}",
                    nearest_gene=nearest_gene, nearest_distance=nearest_dist,
                    n_molqtl_genes=len(molq), molqtl_multiplicity=mult_str, l2g_genes=l2g_str)

    # Tier 3: nearest GENCODE-PC gene
    return dict(assigned_gene=nearest_gene, assignment_tier="nearest",
                assignment_evidence=f"closest|dist={nearest_dist}bp",
                nearest_gene=nearest_gene, nearest_distance=nearest_dist,
                n_molqtl_genes=len(molq), molqtl_multiplicity=mult_str, l2g_genes=l2g_str)


def assign_per_region(df, cache, near):
    """Build the one-gene-per-locus foreground table."""
    near_idx = near.set_index("variant_id")
    rows = []
    for region, g in df.groupby("region"):
        chrom, start, end = region.split("_")
        vids = g["variant_id"].tolist()
        near_sub = near_idx.reindex(vids).reset_index()
        a = assign_one_gene(vids, cache, near_sub)
        rows.append(dict(
            region=region, chrom=chrom, start=int(start), end=int(end),
            n_snps=len(g), mean_bmap=round(pd.to_numeric(g["BMAP"], errors="coerce").mean(), 1),
            **a,
        ))
    cols = ["region", "chrom", "start", "end", "n_snps", "mean_bmap",
            "assigned_gene", "assignment_tier", "assignment_evidence",
            "nearest_gene", "nearest_distance",
            "n_molqtl_genes", "molqtl_multiplicity", "l2g_genes"]
    return pd.DataFrame(rows)[cols].sort_values(["chrom", "start"]).reset_index(drop=True)


# ── Main ───────────────────────────────────────────────────────────────────────

def parse_args():
    """Optional overrides. Defaults reproduce the published 59-region pooled assignment.

    --snp-table lets an alternative TSP set be assigned through identical code; the
    variant cache is shared, so only variants not already seen cost an API call.
    """
    import argparse
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0] if __doc__ else None)
    p.add_argument("--snp-table", default=str(SNP_TABLE))
    p.add_argument("--cache", default=str(CACHE_JSON))
    p.add_argument("--out-variants", default=str(OUT_TSV))
    p.add_argument("--out-regions", default=str(REGION_TSV))
    p.add_argument("--gencode", default=str(GENCODE_BED))
    return p.parse_args()


def main():
    global SNP_TABLE, CACHE_JSON, OUT_TSV, REGION_TSV, GENCODE_BED
    a = parse_args()
    SNP_TABLE, CACHE_JSON = Path(a.snp_table), Path(a.cache)
    OUT_TSV, REGION_TSV, GENCODE_BED = Path(a.out_variants), Path(a.out_regions), Path(a.gencode)

    df = pd.read_csv(SNP_TABLE, sep="\t", dtype=str)
    df["variant_id"] = (
        df["CHROM_HG38"] + "_" + df["POSITION_HG38"] + "_" + df["REF"] + "_" + df["ALT"]
    )
    df = df.drop_duplicates("variant_id").reset_index(drop=True)
    print(f"Variants: {len(df):,}  (regions: {df['region'].nunique()})")

    print("Fetching OpenTargets functional evidence…")
    cache = fetch_variants(df["variant_id"].tolist())

    print("Parsing functional evidence…")
    parsed = pd.DataFrame([parse_variant(v, cache.get(v)) for v in df["variant_id"]])

    print("Tier 3: nearest gene (bedtools closest vs GENCODE v49)…")
    near = nearest_genes(df)

    merged = (
        df[["variant_id", "CHROM_HG38", "POSITION_HG38", "REF", "ALT",
            "BMAP", "Annotation", "Gene_Name", "region"]]
        .merge(parsed, on="variant_id", how="left")
        .merge(near,   on="variant_id", how="left")
    )

    # ── One-gene-per-locus foreground (L2G -> supported molQTL -> nearest) ─────
    region_tbl = assign_per_region(df, cache, near)

    # Per-variant assignment INHERITS its locus decision, so the per-variant and per-locus
    # tables always agree. Deciding per variant instead lets a near-zero-PIP eQTL outrank the
    # aggregated L2G/molQTL evidence, which sends GYPB to the GUSBP5 pseudogene and MTRR to an
    # unnamed ENSG model.
    reg_map = region_tbl.set_index("region")[
        ["assigned_gene", "assignment_tier", "assignment_evidence"]]
    merged = merged.merge(reg_map, on="region", how="left")

    col_order = [
        "variant_id", "rsid", "CHROM_HG38", "POSITION_HG38", "REF", "ALT",
        "BMAP", "region", "Annotation", "most_severe",
        "assigned_gene", "assignment_tier", "assignment_evidence",
        "in_ot_index", "n_credible_sets",
        "best_molqtl_gene", "best_molqtl_pip", "best_molqtl_type", "best_molqtl_tissue",
        "molqtl_genes", "molqtl_details", "gwas_l2g_genes",
        "nearest_gene", "nearest_distance",
        "Gene_Name",  # original snpEff annotation, for comparison
    ]
    merged = merged[col_order]
    merged.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"\nSaved per-variant table (inherits per-locus assignment): {OUT_TSV}")

    region_tbl.to_csv(REGION_TSV, sep="\t", index=False)
    print(f"Saved per-locus foreground: {REGION_TSV}")

    tier_counts = region_tbl["assignment_tier"].value_counts().to_dict()
    print(f"\n  {len(region_tbl)} loci  ->  tiers: {tier_counts}")
    print(f"  unique foreground genes: {region_tbl['assigned_gene'].nunique()}")

    print("\nForeground genes:")
    print("  " + ", ".join(sorted(region_tbl["assigned_gene"].dropna().unique())))


if __name__ == "__main__":
    main()
