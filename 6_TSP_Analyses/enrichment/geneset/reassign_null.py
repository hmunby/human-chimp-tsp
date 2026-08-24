"""
Re-assign an already-sampled, already-fetched null under the CURRENT assignment
rule, WITHOUT re-querying OpenTargets. Reads a read-only snapshot of the null
cache so it is safe to run alongside an in-progress build_background run.

Sampling is deterministic (seed 0), so regenerating the n_rep=200 manifest
reproduces exactly the loci/SNPs that the 200-rep run cached.

Output: null_gene_sets_newlogic.tsv  (one gene-set per replicate, new rule)
"""

import json, shutil, tempfile, time
from pathlib import Path
import pandas as pd

import build_background as bb
import opentargets_assign as ota

N_REP = 200
OUT = bb.GA_DIR / "null_gene_sets_newlogic.tsv"


def snapshot_cache():
    """Copy the live cache and parse it, retrying if caught mid-write."""
    for _ in range(10):
        tmp = Path(tempfile.mkdtemp()) / "cache.json"
        shutil.copy(bb.NULL_CACHE, tmp)
        try:
            return json.loads(tmp.read_text())
        except json.JSONDecodeError:
            time.sleep(1)
    raise RuntimeError("could not get a clean cache snapshot")


def main():
    fg, _ = bb.load_foreground()
    pool, by_chrom = bb.load_pool()
    manifest = bb.sample_null_loci(fg, by_chrom, pool, N_REP, seed=0)
    print(f"regenerated {len(manifest):,} null loci (seed 0, n_rep {N_REP})")

    # explode to SNPs
    recs = []
    for row in manifest.itertuples():
        uid = f"{row.fg_region}|rep{row.replicate}"
        for p in str(row.null_positions).split(","):
            recs.append((uid, row.replicate, row.null_chrom, int(p)))
    snps = pd.DataFrame(recs, columns=["locus_uid", "replicate", "chrom", "pos"])

    alleles = bb.fetch_alleles(snps[["chrom", "pos"]].drop_duplicates())
    def vid(r):
        a = alleles.get((r.chrom, r.pos))
        return f"{r.chrom.replace('chr','')}_{r.pos}_{a[0]}_{a[1]}" if a else None
    snps["variant_id"] = [vid(r) for r in snps.itertuples()]
    snps = snps.dropna(subset=["variant_id"])

    cache = snapshot_cache()
    hit = snps["variant_id"].isin(cache).mean()
    print(f"cache snapshot: {len(cache):,} SNPs | 200-rep variants covered: {hit*100:.1f}%")

    snp_df = snps.drop_duplicates("variant_id").copy()
    snp_df["CHROM_HG38"] = snp_df["chrom"].str.replace("chr", "", regex=False)
    snp_df["POSITION_HG38"] = snp_df["pos"]
    near = ota.nearest_genes(snp_df[["CHROM_HG38", "POSITION_HG38", "variant_id"]])
    near_idx = near.set_index("variant_id")

    rows = []
    for uid, g in snps.groupby("locus_uid"):
        vids = g["variant_id"].tolist()
        a = ota.assign_one_gene(vids, cache, near_idx.reindex(vids).reset_index())
        rows.append(dict(replicate=int(g["replicate"].iloc[0]),
                         assigned_gene=a["assigned_gene"], tier=a["assignment_tier"]))
    nl = pd.DataFrame(rows)

    sets = (nl.groupby("replicate")["assigned_gene"]
              .apply(lambda s: ",".join(sorted(set(s.dropna())))).reset_index())
    sets.to_csv(OUT, sep="\t", index=False)
    print(f"\nsaved: {OUT}")
    print(f"null tiers (new rule): {nl['tier'].value_counts().to_dict()}")


if __name__ == "__main__":
    main()
