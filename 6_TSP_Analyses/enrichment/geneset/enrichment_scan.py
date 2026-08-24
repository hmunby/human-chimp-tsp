"""
Empirical-null enrichment scan across GO / KEGG / Reactome / Hallmark categories.

For every category that contains >= MIN_K of the foreground genes, test enrichment
TWO ways and compare:

  unmatched : foreground vs all autosomal PC genes (hypergeometric) — what a default
              whole-genome DAVID background gives.
  matched   : foreground vs the BMAP+structure+rule matched null gene-sets (empirical).

A category enriched vs the unmatched null but NOT vs the matched null is an
ASCERTAINMENT ARTIFACT (like glycoprotein). A category enriched vs BOTH is a
candidate real signal that survives ascertainment matching.

Outputs:
  enrichment_scan_results.tsv  - every tested category, both tests, BH-FDR, verdict
"""

from pathlib import Path
import argparse
import glob
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent

# The foreground gene set produced by 4_TSPRegions/3_GeneAssignment/.
REGION_TSV    = HERE / "../../../4_SINGER_GenealogicalReconstruction_TSP/4_TSPRegions/" \
                       "3_GeneAssignment/region_gene_assignments.pooled.tsv"
NULL_SETS_TSV = HERE / "null_gene_sets.tsv"
GENESET_DIR   = HERE / "genesets"
# GENCODE v49 protein-coding autosomal genes as BED with gene names; not bundled.
GENCODE_BED   = HERE / "resources/gencode.v49.annotation.genes.protein_coding.autosomes." \
                       "gene_names.bed"
OUT_TSV = HERE / "enrichment_scan_results.tsv"

MIN_K = 2     # default; override with --min-k. k=2 turns the foreground into a
              # noise generator (any co-occurring PAIR gives a large fold), so k>=3 is
              # the sensible primary now that the test has power at 2,000 replicates.


def load_gmt(path):
    """library_name -> {category: set(genes)}"""
    cats = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            name, _desc, *genes = parts
            cats[name] = {g.split(",")[0].strip().upper() for g in genes if g.strip()}
    return cats


def bh_fdr(pvals):
    p = np.asarray(pvals, float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--null", default=str(NULL_SETS_TSV),
                    help="per-replicate null gene-sets TSV")
    ap.add_argument("--regions", default=str(REGION_TSV),
                    help="foreground region_gene_assignments.tsv; override to test an "
                         "alternative TSP set through identical code")
    ap.add_argument("--gencode", default=str(GENCODE_BED))
    ap.add_argument("--min-k", type=int, default=MIN_K,
                    help="minimum foreground genes for a category to be tested")
    ap.add_argument("--out", default=str(OUT_TSV))
    args = ap.parse_args()
    globals()["GENCODE_BED"] = args.gencode

    fg = [g.upper() for g in pd.read_csv(args.regions, sep="\t")["assigned_gene"].dropna()]
    fg_set = set(fg)
    # Number of DRAWS for the hypergeometric must be the number of DISTINCT genes tested, since
    # k_obs is computed from fg_set. Using one draw per region instead would over-count the genes
    # that are assigned to more than one region, inflating the expectation and making
    # fold_unmatched/p_unmatched conservative.
    n = len(fg_set)

    universe = {g.upper() for g in pd.read_csv(
        GENCODE_BED, sep="\t", header=None, names=["c", "s", "e", "g"])["g"]}
    N = len(universe)

    sets = pd.read_csv(args.null, sep="\t")
    print(f"null gene-sets: {args.null}")
    null_sets = [set(str(r["assigned_gene"]).upper().split(",")) for _, r in sets.iterrows()]
    n_rep = len(null_sets)
    null_sz = np.array([len(x) for x in null_sets])   # null set sizes, for the fraction statistic
    print(f"Foreground: {n} genes | universe: {N:,} | matched null replicates: {n_rep}")

    rows = []
    for gmt in sorted(glob.glob(str(GENESET_DIR / "*.gmt"))):
        lib = Path(gmt).stem
        cats = load_gmt(gmt)
        tested = 0
        for cat, genes in cats.items():
            k_obs = len(fg_set & genes)
            if k_obs < args.min_k:
                continue
            tested += 1
            # unmatched: hypergeometric vs all PC genes
            K = len(genes & universe)
            p_unmatched = float(stats.hypergeom.sf(k_obs - 1, N, K, n)) if K else 1.0
            exp_unmatched = n * K / N
            # matched: empirical vs null gene-sets
            null_k = np.array([len(s & genes) for s in null_sets])
            p_matched = float((1 + (null_k >= k_obs).sum()) / (1 + n_rep))
            exp_matched = null_k.mean()
            # FRACTION statistic. The count statistic above compares k_obs against null COUNTS,
            # which is only valid if the null gene sets are the same size as the foreground. Here
            # they are close (mean 48.6 vs 46, ratio 1.057) so the bias is ~6%, unlike the
            # uncoalesced analysis where the ratio was 1.30-1.49 and the count statistic reversed
            # the conclusion. Reported alongside so the two are comparable.
            null_frac = null_k / null_sz
            obs_frac = k_obs / n
            fold_matched_frac = obs_frac / null_frac.mean() if null_frac.mean() else np.inf
            p_matched_frac = float((1 + (null_frac >= obs_frac).sum()) / (1 + n_rep))
            rows.append(dict(
                library=lib, category=cat, k_obs=k_obs,
                fg_genes=",".join(sorted(fg_set & genes)),
                exp_unmatched=round(exp_unmatched, 2),
                fold_unmatched=round(k_obs / exp_unmatched, 2) if exp_unmatched else np.inf,
                p_unmatched=p_unmatched,
                exp_matched=round(exp_matched, 2),
                fold_matched=round(k_obs / exp_matched, 2) if exp_matched else np.inf,
                p_matched=p_matched,
                fold_matched_frac=round(fold_matched_frac, 3),
                p_matched_frac=p_matched_frac,
            ))
        print(f"  {lib}: {tested} categories with >={args.min_k} foreground genes")

    res = pd.DataFrame(rows)
    if res.empty:
        print("No categories met the threshold.")
        return

    # BH-FDR across all tested categories (pooled), for each test
    res["p_unmatched_BH"] = bh_fdr(res["p_unmatched"])
    res["p_matched_BH"]   = bh_fdr(res["p_matched"])
    res["p_matched_frac_BH"] = bh_fdr(res["p_matched_frac"])

    # The verdict uses the FRACTION statistic and its BH value. The count statistic compares
    # k_obs against null COUNTS without accounting for null gene-set size (mean 48.6 vs 46
    # observed); a held-out calibration on 2,000 replicates puts the fraction version at 4.48%
    # rejection against a nominal 5% and the count version at 3.49%, i.e. both conservative but
    # the count one needlessly so. The count columns are retained for comparison.
    def verdict(r):
        sig_un = r["p_unmatched_BH"] < 0.05 and r["fold_unmatched"] > 1
        sig_ma = r["p_matched_frac_BH"] < 0.05 and r["fold_matched_frac"] > 1
        if sig_ma:
            return "real_survives_matching"
        if sig_un and not sig_ma:
            return "ascertainment_artifact"
        return "ns"
    res["verdict"] = res.apply(verdict, axis=1)

    res = res.sort_values(["p_matched_frac", "p_unmatched"]).reset_index(drop=True)
    res.to_csv(args.out, sep="\t", index=False)

    print(f"\nSaved: {args.out}  ({len(res)} tested categories)")
    print("\nverdict counts:", res["verdict"].value_counts().to_dict())

    surv = res[res["verdict"] == "real_survives_matching"].copy()
    print(f"\n── Categories enriched vs the MATCHED null (fraction statistic, BH<0.05): {len(surv)} ──")
    cols = ["library", "category", "k_obs", "fold_matched_frac", "p_matched_frac",
            "p_matched_frac_BH", "fg_genes"]
    with pd.option_context("display.width", 260, "display.max_colwidth", 60):
        print(surv[cols].to_string(index=False))


if __name__ == "__main__":
    main()
