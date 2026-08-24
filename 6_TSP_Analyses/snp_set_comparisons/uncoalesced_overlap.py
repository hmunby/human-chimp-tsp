#!/usr/bin/env python3
"""Overlap of each SNP set with regions uncoalesced in both species (Figure 2).

For the four nested SNP sets (1000GP MAF >= 5% > all shared > old shared >= 4 Mya > TSP),
the percentage of SNPs falling in the reciprocal (human x chimp) uncoalesced segments from
stage 5, against a size-matched null.

Null: the observed count for each set is compared to `N_PERM` random subsamples of the same
size drawn from the 1000GP MAF >= 5% pool. This matches set size and the common-SNP
ascertainment only -- it is deliberately NOT matched on BMAP, age or mutation rate, so the
fold-enrichments are the raw excess over common variation. p = (#null >= obs + 1) / (N + 1),
so the floor with 1000 permutations is 9.99e-4. The 1000GP row is the null pool itself and
so has no expectation or p-value.

Inputs:  the reciprocal intersect from ../../5_UncoalescedRegions/, the 1000GP and
         all-shared hg38 SNP beds (large, not bundled), and the two bundled sets in misc_files/.
Outputs: snp_set_uncoal_overlap.{txt,png,no_expected.png}
"""
import os
import subprocess

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))

# ---- inputs (edit) -------------------------------------------------------------------------
# Reciprocal uncoalesced segments (stage 5, 2_UncoalescedSegments)
UNCOAL = os.path.join(HERE, "resources/expected_uncoal_intersect.whole_genome.exact.bed")

# Large SNP beds (hg38, CHROM/START/END), not bundled -- see README
TGP_BED = os.path.join(HERE, "resources/tgp_snps_hg38.bed")                  # 153 MB
SHARED_BED = os.path.join(HERE, "resources/shared_snps_chrom_pos_hg38.sorted.bed")   # 1.6 MB

# Bundled sets. The old-shared bed is the 860 candidate set lifted to hg38; 858 of the 860 lift
# (9:136143120 and 14:106529056 do not).
OLD_SHARED_BED = os.path.join(HERE, "misc_files/old_shared_snps_860.hg38.sorted.bed")
TSP_BED = os.path.join(HERE, "misc_files/tsp_snps_pooled.hg38.bed")

BEDTOOLS = "bedtools"
N_PERM = 1000
SEED = 0
# --------------------------------------------------------------------------------------------

SETS = [
    ("1000GP (MAF ≥ 5%)", TGP_BED, "#7F7F7F"),
    ("All shared", SHARED_BED, "#1F77B4"),
    ("Old shared (≥4 Mya)", OLD_SHARED_BED, "#9467BD"),
    ("TSP", TSP_BED, "#FF0000"),
]


def n_overlap(df):
    """Number of SNPs in `df` overlapping the uncoalesced segments."""
    p = subprocess.run([BEDTOOLS, "intersect", "-u", "-a", "-", "-b", UNCOAL],
                       input=df.to_csv(sep="\t", header=False, index=False),
                       capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(p.stderr)
    return 0 if not p.stdout else p.stdout.count("\n")


def null_distribution(tgp_df, n, n_perm=None, rng=None):
    """`n_perm` overlap counts for random size-`n` subsamples of the 1000GP pool."""
    n_perm = N_PERM if n_perm is None else n_perm
    rng = np.random.default_rng(SEED) if rng is None else rng
    null = np.empty(n_perm)
    for i in range(n_perm):
        samp = tgp_df.sample(n=n, replace=False, random_state=int(rng.integers(2**31 - 1)))
        null[i] = n_overlap(samp)
    return null


def overlap_row(name, bed, color, tgp_df, rng, is_pool=False):
    """One SNP set's observed overlap and its size-matched null."""
    df = pd.read_csv(bed, sep="\t", header=None, names=["c", "s", "e"])
    N = len(df)
    obs = n_overlap(df)
    if is_pool:
        exp_mean = exp_sd = float(obs)      # the null pool itself
        pval = np.nan
    else:
        null = null_distribution(tgp_df, N, rng=rng)
        exp_mean, exp_sd = null.mean(), null.std()
        pval = (np.sum(null >= obs) + 1) / (N_PERM + 1)
    return dict(snp_set=name, N=N, observed=obs, obs_pct=100 * obs / N,
                exp_mean=exp_mean, exp_pct=100 * exp_mean / N,
                exp_ci=100 * 1.96 * exp_sd / N,
                fold=(obs / exp_mean) if exp_mean > 0 else np.nan, p=pval, color=color)


def compute(sets=None):
    """The four-row result table. `sets` overrides SETS, e.g. to swap the TSP bed."""
    rng = np.random.default_rng(SEED)
    tgp_df = pd.read_csv(TGP_BED, sep="\t", header=None, names=["c", "s", "e"])
    rows = [overlap_row(name, bed, color, tgp_df, rng, is_pool=name.startswith("1000GP"))
            for name, bed, color in (sets or SETS)]
    return pd.DataFrame(rows)


def bars(ax, res):
    x = np.arange(len(res))
    width = 0.6
    ax.bar(x, res["obs_pct"], color=res["color"], width=width, zorder=2)
    for xi, r in enumerate(res.itertuples()):
        ax.text(xi, r.obs_pct + max(res["obs_pct"]) * 0.02, f"{r.obs_pct:.1f}%",
                ha="center", va="bottom", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels(res["snp_set"], rotation=20, ha="right")
    ax.set_ylabel("Percentage of SNPs in segments\nuncoalesced in both species", fontweight="bold")
    ax.set_ylim(0, max(res["obs_pct"]) * 1.18)
    ax.spines[["top", "right"]].set_visible(False)
    return x, width


def plot(res, out_png, out_png_no_expected):
    # with the expected-value overlay
    fig, ax = plt.subplots(figsize=(5.4, 4.4), dpi=300)
    x, width = bars(ax, res)
    for xi, r in enumerate(res.itertuples()):
        if not r.snp_set.startswith("1000GP"):
            ax.hlines(r.exp_pct, xi - width / 2, xi + width / 2, color="black", lw=2.5, zorder=5)
            ax.errorbar(xi, r.exp_pct, yerr=r.exp_ci, fmt="none", ecolor="black",
                        elinewidth=1.2, capsize=3, zorder=5)
    leg = ax.legend(handles=[Line2D([0], [0], color="black", lw=2.5,
                                    label="Expected (1000GP subsample, 95% CI)")],
                    fontsize=8, loc="upper left", frameon=True)
    leg.get_frame().set_edgecolor("black")
    leg.get_frame().set_linewidth(0.8)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300)

    # observed only
    fig2, ax2 = plt.subplots(figsize=(5.4, 4.4), dpi=300)
    bars(ax2, res)
    fig2.tight_layout()
    fig2.savefig(out_png_no_expected, dpi=300)


if __name__ == "__main__":
    import argparse
    _ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    _ap.add_argument("--tsp-bed", default=TSP_BED,
                     help="alternative TSP set (hg38 BED); defaults to the reported pooled set")
    _ap.add_argument("--out-prefix", default=os.path.join(HERE, "snp_set_uncoal_overlap"))
    _a = _ap.parse_args()

    _sets = [(n, _a.tsp_bed if n == "TSP" else b, c) for n, b, c in SETS]
    res = compute(_sets)
    res.drop(columns="color").to_csv(f"{_a.out_prefix}.txt", sep="\t", index=False)
    plot(res, f"{_a.out_prefix}.png", f"{_a.out_prefix}.no_expected.png")
    print(res.drop(columns="color").to_string(index=False))
