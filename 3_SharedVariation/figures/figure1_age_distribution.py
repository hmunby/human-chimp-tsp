#!/usr/bin/env python3
"""Figure 1 - age distribution of shared SNPs and mutation rate vs age.

Panel A: percentage of SNPs older than a given age, for shared human-chimp SNPs vs all common
         1000GP (MAF>=5%) SNPs, with the human-chimp split shaded at 5-6 Mya.
Panel B: mean Roulette mutation rate of shared SNPs as a function of a minimum YRI age, vs the
         genome-wide 1000GP (MAF>=5%) mean.

Ages are the YRI estimates throughout. The paper figure is MHC-excluded; both versions are written:
    figures/figure1_age_distribution.png            (MHC excluded; the paper figure)
    figures/figure1_age_distribution.with_mhc.png   (MHC included; alternate)

Run:  python figure1_age_distribution.py
Paths and constants are in config.py; the shaded-band / log-tick helpers are in plotting.py.

Run with `common/` on PYTHONPATH (this imports `config` and `plotting`):

    PYTHONPATH=../../common python figure1_age_distribution.py

Input paths and constants come from common/config.py. Large raw inputs are not bundled; the
small derived tables the script reads are.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import config
from plotting import gradient_vspan, plain_log_fmt


def load_panel_a(tgp_path, shared_path):
    """Expected-SNPs-by-age curves; return (shared_df, tgp_df, shared_total, tgp_total)."""
    tgp = pd.read_csv(tgp_path, sep="\t")
    shared = pd.read_csv(shared_path, sep="\t")
    tgp_total = tgp.loc[tgp["age"] == 0.0, "expected_snps"].values[0]
    shared_total = shared.loc[shared["age"] == 0.0, "expected_snps"].values[0]
    for d, tot in [(tgp, tgp_total), (shared, shared_total)]:
        d["fraction"] = d["expected_snps"] / tot
        d["age_mya"] = d["age"] / 1e6      # the 'age' column is already in years
    return shared, tgp, shared_total, tgp_total


def compute_panel_b(ann_path, exclude_mhc):
    """Mean shared-SNP mutation rate as a function of a minimum YRI mid-age; return (ages, mean, ci).

    Binned by MID_AGE_YRI, the YRI Relate age estimates used throughout.
    """
    df = pd.read_csv(ann_path, sep="\t", low_memory=False)
    if exclude_mhc:
        c, s, e = config.MHC_B37
        df = df[~((df["CHROM"].astype(str) == c) & df["POS"].between(s, e))]
    # the annotation table writes "." for SNPs Relate could not date
    for col in ("MR", "LOWER_AGE_YRI", "UPPER_AGE_YRI"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["MID_AGE_YRI"] = df[["UPPER_AGE_YRI", "LOWER_AGE_YRI"]].mean(axis=1)
    ages = list(range(0, 9))
    mean, ci = [], []
    for age in ages:
        sub = df[df["MID_AGE_YRI"] >= age * 1e6 / config.GEN_YEARS]
        mean.append(sub["MR"].mean())
        ci.append(1.96 * sub["MR"].std() / np.sqrt(len(sub)))
    return ages, mean, ci


def make_figure(exclude_mhc, out_png):
    if exclude_mhc:
        tgp_p, shared_p = config.FIG1_TGP_EXPECTED_NOMHC, config.FIG1_SHARED_EXPECTED_NOMHC
    else:
        tgp_p, shared_p = config.FIG1_TGP_EXPECTED, config.FIG1_SHARED_EXPECTED
    tag = "MHC excluded" if exclude_mhc else "MHC included"

    shared_df, tgp_df, shared_total, tgp_total = load_panel_a(tgp_p, shared_p)
    ages, mr_mean, mr_ci = compute_panel_b(config.FIG1_SHARED_ANN, exclude_mhc)

    print(f"[{tag}]")
    for cutoff in [4, 5]:
        row = shared_df.loc[shared_df["age_mya"] >= cutoff]
        n, frac = row["expected_snps"].values[0], row["fraction"].values[0]
        print(f"  Expected shared SNPs older than {cutoff} Mya: {n:.1f} ({frac*100:.3f}% of total)")

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(10, 4), dpi=400)

    # --- Panel A ---
    ax_a.plot(shared_df["age_mya"], shared_df["fraction"] * 100,
              label=f"Shared SNPs (n={int(shared_total):,})", color="#1f77b4")
    ax_a.plot(tgp_df["age_mya"], tgp_df["fraction"] * 100,
              label=f"1000GP (MAF ≥ 5%) (n={int(tgp_total):,})", color="red")
    ax_a.set_xlabel("Age (millions of years)", fontweight="bold")
    ax_a.set_ylabel("Percentage of SNPs (%)", fontweight="bold")
    ax_a.set_yscale("log")
    ax_a.set_xlim(0, 10)
    ax_a.set_ylim(1e-2, 100)
    ax_a.yaxis.set_major_formatter(plain_log_fmt)
    hc_patch = gradient_vspan(ax_a, core=config.HC_SPLIT_MYA, fade=0.5)
    h, l = ax_a.get_legend_handles_labels()
    h.append(hc_patch); l.append("Human-Chimp split")
    ax_a.legend(h, l, loc="lower left", frameon=False)
    ax_a.text(-0.12, 1.05, "A", transform=ax_a.transAxes, fontsize=14, fontweight="bold")

    # --- Panel B ---
    ax_b.errorbar(ages, mr_mean, yerr=mr_ci, label="Shared SNPs", marker="o", elinewidth=2)
    ax_b.axhline(y=config.TGP_MEAN_MR, color="r", linestyle="--", label="All 1000GP (MAF ≥ 5%)")
    ax_b.set_xlabel("Age (millions of years)", fontweight="bold")
    ax_b.set_ylabel("Mean mutation rate", fontweight="bold")
    ax_b.set_ylim(0, 1.2)
    ax_b.set_xlim(0, 8.1)
    ax_b.legend(frameon=False)
    ax_b.text(-0.12, 1.05, "B", transform=ax_b.transAxes, fontsize=14, fontweight="bold")

    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)
    print(f"  wrote {out_png}")


if __name__ == "__main__":
    os.makedirs(config.FIG_OUTDIR, exist_ok=True)
    # canonical paper figure = MHC-excluded (both panels); MHC-included version kept as alternate
    make_figure(exclude_mhc=True,
                out_png=os.path.join(config.FIG_OUTDIR, "figure1_age_distribution.png"))
    make_figure(exclude_mhc=False,
                out_png=os.path.join(config.FIG_OUTDIR, "figure1_age_distribution.with_mhc.png"))
