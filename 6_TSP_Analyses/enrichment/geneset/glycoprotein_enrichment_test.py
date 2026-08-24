"""
Glycoprotein (UniProt KW-0325) enrichment in the test loci, tested two ways:

  (1) NAIVE / UNMATCHED null  - foreground vs all autosomal protein-coding genes
      (what a default whole-genome DAVID background does). Hypergeometric +
      empirical random-gene draws.

  (2) MATCHED null            - foreground vs BMAP + locus-structure + assignment-
      rule matched null loci (the empirical-null gene sets from build_background.py).

The contrast is the headline result: the apparent glycoprotein enrichment relative
to all genes is explained by ascertainment (test loci sit in high-BMAP, gene-sparse
regions which are intrinsically glycoprotein-rich); it disappears against the
matched null.

Outputs:
  glycoprotein_enrichment_results.tsv   - the numbers
  glycoprotein_enrichment_results.md    - human-readable summary
"""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent

REGION_TSV    = HERE / "../../../4_SINGER_GenealogicalReconstruction_TSP/4_TSPRegions/" \
                       "3_GeneAssignment/region_gene_assignments.pooled.tsv"
NULL_SETS_TSV = HERE / "null_gene_sets.tsv"
# GENCODE v49 protein-coding autosomal genes as BED with gene names; not bundled.
GENCODE_BED   = HERE / "resources/gencode.v49.annotation.genes.protein_coding.autosomes." \
                       "gene_names.bed"
GLYCO_TSV = HERE / "../glycoprotein/human_glycoproteins_KW0325_with_coords.tsv"
OUT_TSV = HERE / "glycoprotein_enrichment_results.tsv"
OUT_MD  = HERE / "glycoprotein_enrichment_results.md"

N_DRAWS = 100_000
SEED    = 0



def _cli():
    """Optional overrides. Defaults reproduce the published result.

    --regions points the test at an alternative TSP set's gene assignment, and --null at
    a null resampled to that set's size; both are used by the Ne-sensitivity comparison.
    """
    import argparse
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1] if __doc__ else None)
    p.add_argument("--regions", default=str(REGION_TSV))
    p.add_argument("--null", default=str(NULL_SETS_TSV))
    p.add_argument("--gencode", default=str(GENCODE_BED))
    p.add_argument("--glyco", default=str(GLYCO_TSV))
    p.add_argument("--out-tsv", default=str(OUT_TSV))
    p.add_argument("--out-md", default=str(OUT_MD))
    a = p.parse_args()
    g = globals()
    g["REGION_TSV"], g["NULL_SETS_TSV"] = Path(a.regions), Path(a.null)
    g["GENCODE_BED"], g["GLYCO_TSV"] = Path(a.gencode), Path(a.glyco)
    g["OUT_TSV"], g["OUT_MD"] = Path(a.out_tsv), Path(a.out_md)


def main():
    glyco = set(pd.read_csv(GLYCO_TSV, sep="\t")["hgnc_symbol"].dropna())
    universe = list(pd.read_csv(GENCODE_BED, sep="\t", header=None,
                                names=["c", "s", "e", "g"])["g"])
    N, K = len(universe), len(set(universe) & glyco)
    bg_rate = K / N

    fg = pd.read_csv(REGION_TSV, sep="\t")
    fg_genes = sorted(set(fg["assigned_gene"].dropna()))   # DISTINCT genes: 46, not 49 rows.
    # The null sets are deduplicated, so counting the foreground with duplicates (IGFBP7, MUC7
    # each appear twice) would inflate the observed count relative to the null.
    n = len(fg_genes)
    obs = sum(1 for g in fg_genes if g in glyco)
    fg_gly = [g for g in fg_genes if g in glyco]

    rng = np.random.default_rng(SEED)

    # ── (1) unmatched / all-genes null ────────────────────────────────────────
    p_hyper = float(stats.hypergeom.sf(obs - 1, N, K, n))
    draws = np.array([
        sum(1 for g in rng.choice(universe, n, replace=False) if g in glyco)
        for _ in range(N_DRAWS)
    ])
    p_unmatched_emp = float((draws >= obs).mean())
    fold_unmatched = (obs / n) / bg_rate

    # ── (2) matched null (empirical) ──────────────────────────────────────────
    sets = pd.read_csv(NULL_SETS_TSV, sep="\t")
    null_k = np.array([
        sum(1 for g in str(row["assigned_gene"]).split(",") if g in glyco)
        for _, row in sets.iterrows()
    ])
    n_rep = len(null_k)
    null_sz = np.array([len([g for g in str(row["assigned_gene"]).split(",") if g])
                        for _, row in sets.iterrows()])
    matched_mean = null_k.mean()
    # COUNT statistic (original): obs vs null COUNTS. Valid only if the null sets are the same size
    # as the foreground; they average 48.6 against 46, so it is mildly biased against enrichment.
    p_matched = float((1 + (null_k >= obs).sum()) / (1 + n_rep))   # +1 smoothing
    fold_matched = (obs / n) / (matched_mean / n) if matched_mean > 0 else np.nan
    # FRACTION statistic: compares composition, which is what the test is asking. A held-out
    # calibration on 2,000 replicates puts this at 4.48% rejection against a nominal 5% (the count
    # version 3.49%), i.e. both conservative, the fraction version less needlessly so.
    null_frac = null_k / null_sz
    obs_frac_v = obs / n
    fold_matched_frac = obs_frac_v / null_frac.mean() if null_frac.mean() > 0 else np.nan
    p_matched_frac = float((1 + (null_frac >= obs_frac_v).sum()) / (1 + n_rep))

    # ── results table ─────────────────────────────────────────────────────────
    rows = [
        dict(test="unmatched_hypergeometric", null_source="all autosomal PC genes",
             n_foreground=n, obs_glyco=obs, obs_frac=round(obs / n, 3),
             null_frac=round(bg_rate, 3), fold=round(fold_unmatched, 2),
             p_value=round(p_hyper, 4), n_null=N),
        dict(test="unmatched_empirical", null_source=f"{N_DRAWS} random {n}-gene draws",
             n_foreground=n, obs_glyco=obs, obs_frac=round(obs / n, 3),
             null_frac=round(draws.mean() / n, 3), fold=round(fold_unmatched, 2),
             p_value=round(p_unmatched_emp, 4), n_null=N_DRAWS),
        dict(test="matched_empirical_COUNT",
             null_source="BMAP+structure+rule matched null loci (count statistic)",
             n_foreground=n, obs_glyco=obs, obs_frac=round(obs / n, 3),
             null_frac=round(matched_mean / n, 3), fold=round(fold_matched, 2),
             p_value=round(p_matched, 4), n_null=n_rep),
        dict(test="matched_empirical_FRACTION",
             null_source="BMAP+structure+rule matched null loci (fraction statistic)",
             n_foreground=n, obs_glyco=obs, obs_frac=round(obs_frac_v, 3),
             null_frac=round(float(null_frac.mean()), 3), fold=round(fold_matched_frac, 2),
             p_value=round(p_matched_frac, 4), n_null=n_rep),
    ]
    res = pd.DataFrame(rows)
    res.to_csv(OUT_TSV, sep="\t", index=False)

    # ── markdown summary ──────────────────────────────────────────────────────
    md = f"""# Glycoprotein (UniProt KW-0325) enrichment in the test loci

**Foreground:** {n} loci -> {n} genes, of which **{obs} are glycoproteins ({100*obs/n:.0f}%)**.
Glycoprotein foreground genes: {', '.join(fg_gly)}

**Genome-wide rate:** {K:,}/{N:,} autosomal protein-coding genes are glycoproteins ({100*bg_rate:.1f}%).

| Null model | Null glyco % | Fold | p |
|---|---|---|---|
| Unmatched — all PC genes (hypergeometric) | {100*bg_rate:.0f}% | {fold_unmatched:.2f}x | {p_hyper:.4f} |
| Unmatched — {N_DRAWS:,} random {n}-gene draws | {100*draws.mean()/n:.0f}% | {fold_unmatched:.2f}x | {p_unmatched_emp:.4f} |
| Matched — count statistic ({n_rep} reps) | {100*matched_mean/n:.0f}% | {fold_matched:.2f}x | {p_matched:.4f} |
| **Matched — fraction statistic ({n_rep} reps)** | **{100*null_frac.mean():.0f}%** | **{fold_matched_frac:.2f}x** | **{p_matched_frac:.4f}** |

## Interpretation
Against an unmatched whole-genome background the loci look enriched
({100*obs/n:.0f}% vs {100*bg_rate:.0f}%, {fold_unmatched:.2f}x). The matched null itself runs at
{100*null_frac.mean():.0f}% — the loci are ascertained in high-BMAP, gene-sparse regions that are
intrinsically ~{null_frac.mean()/bg_rate:.1f}x richer. Against that matched null the fold is
{fold_matched_frac:.2f}x at p = {p_matched_frac:.4f} ({n_rep} replicates), so read the matched row,
not the unmatched one.

Counting note: the foreground is {n} DISTINCT genes, not the region rows,
double-counting genes assigned to two regions.

Generated by glycoprotein_enrichment_test.py.
"""
    OUT_MD.write_text(md)

    print(res.to_string(index=False))
    print(f"\nSaved: {OUT_TSV}")
    print(f"Saved: {OUT_MD}")


if __name__ == "__main__":
    _cli()
    main()
