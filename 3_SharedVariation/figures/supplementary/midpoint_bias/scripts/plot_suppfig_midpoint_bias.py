#!/usr/bin/env python3
"""Supplementary Figure 3 - bias in the mid-branch (midpoint) age estimator.

Regenerates the figure `SuppFig_midpoint_bias_generations_100k{,.notitle}.png`:
a scatter of the midpoint age estimate vs the true mutation age (both in generations)
for a random subsample of mutations, the y = x line, and a loess of the midpoint
estimate on the true age computed over all simulated mutations.

Underlying simulation (scripts/simulate_neutral.py + Snakefile), verified as the source:
  msprime.sim_ancestry(samples=2000)  -> 2000 diploid individuals = 4000 haploid genomes
  population_size Ne = 10,000 ; mu = 1e-8 /site/gen ; 1 Mb per replicate, no recombination
  mutations kept if derived-allele frequency >= 1% and not on the root branch
  ~50,000 replicates -> 100,090,793 mutations (results/neutral_sim_combined.pkl)

Two modes:
  (default) plot from the archived cache  scripts/midpoint_curve_cache.npz
            -> reproduces the published figure exactly (8000-point subsample + stored loess).
  --rebuild recompute the 8000-point subsample and the loess from the raw combined pickle
            (results/neutral_sim_combined.pkl) and refresh the cache, then plot.
            APPROXIMATE: the archived figure's exact loess parameters were not saved, so the
            rebuilt curve is close but not identical (crossover ~2.3-3.2 Ne vs the archived
            3.14 Ne, depending on LOESS_FRAC). The loess here is a lowess-smoothed log-binned
            conditional mean E[midpoint | true_age] -- a plain lowess on the raw points fails
            badly because true_age is very young-skewed (crossover collapses to <1 Ne). Use the
            default (from-cache) mode to reproduce the published figure exactly.

Usage:
  python plot_suppfig_midpoint_bias.py                # plot from cache
  python plot_suppfig_midpoint_bias.py --rebuild      # recompute from raw data, then plot
"""
import argparse, os, pickle
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "midpoint_curve_cache.npz")
COMBINED = os.path.join(os.path.dirname(HERE), "results", "neutral_sim_combined.pkl")

# --- figure style (matches the published SuppFig) ---
N_SUB = 8000          # scatter subsample size
AXIS_MAX = 100_000    # generations; the "_100k" axis cap
SEED = 0
SCATTER_C = "#4C72B0"     # steelblue points
LOESS_C = "#C44E52"       # muted red loess
LOESS_FRAC = 0.10         # lowess span for smoothing the binned conditional means (rebuild only)
LOESS_NPTS = 250          # loess evaluation points (rebuild only)


def rebuild_cache():
    """Recompute the 8000-point subsample and the loess from the raw combined pickle."""
    from statsmodels.nonparametric.smoothers_lowess import lowess
    print(f"loading {COMBINED} ...")
    with open(COMBINED, "rb") as f:
        data = pickle.load(f)
    arr = np.asarray(data, dtype=np.int64)          # (N, 2): (true_age, midpoint_age), generations
    true_age, mid_age = arr[:, 0], arr[:, 1]
    print(f"  {len(true_age):,} mutations")

    rng = np.random.default_rng(SEED)
    idx = rng.choice(len(true_age), size=N_SUB, replace=False)
    ax_s, mx_s = true_age[idx], mid_age[idx]

    # The "loess" curve is the conditional mean E[midpoint | true_age] over ALL mutations.
    # NOTE: a plain lowess on the raw points does NOT work here - the true-age distribution is
    # very young-skewed, so a fraction-of-points window spans a huge x-range at old ages and
    # drags the fit down (crossover collapses to <1 Ne). Instead we bin by true age on a log
    # scale, take the mean midpoint per bin (the conditional mean), then lowess-smooth those
    # bin means into a curve. This reproduces the archived crossover (~3.1-3.2 Ne). Verified
    # against the alternatives in the method experiment (binned-mean = 3.22 Ne; raw lowess = 0.6 Ne).
    xmax = float(true_age.max())
    edges = np.logspace(np.log10(max(true_age.min(), 1)), np.log10(xmax), 120)
    ctr = np.sqrt(edges[:-1] * edges[1:])
    bidx = np.digitize(true_age, edges)
    bmean = np.full(len(ctr), np.nan)
    for i in range(1, len(edges)):
        m = bidx == i
        if m.sum() >= 20:
            bmean[i - 1] = mid_age[m].mean()
    ok = ~np.isnan(bmean)
    xvals = np.linspace(true_age.min(), xmax, LOESS_NPTS)
    print(f"  loess = lowess-smoothed log-binned conditional mean ({ok.sum()} bins) "
          f"evaluated at {LOESS_NPTS} grid points ...")
    ly = lowess(bmean[ok], ctr[ok], frac=LOESS_FRAC, xvals=xvals)
    lx = xvals
    np.savez(CACHE, ax_s=ax_s, mx_s=mx_s, lx=lx, ly=ly, xmax=xmax)
    print(f"  wrote {CACHE}")
    return ax_s, mx_s, lx, ly, xmax


def load_cache():
    z = np.load(CACHE)
    return z["ax_s"], z["mx_s"], z["lx"], z["ly"], float(z["xmax"])


def make_figure(ax_s, mx_s, lx, ly, out_prefix):
    # scatter transparency/size scale with the subsample so denser plots stay readable
    n = len(ax_s)
    alpha = float(np.clip(2500.0 / n, 0.05, 0.30))
    msize = 8 if n <= 10_000 else 4
    for titled in (True, False):
        fig, ax = plt.subplots(figsize=(8, 7.5), dpi=300)
        ax.scatter(ax_s, mx_s, s=msize, alpha=alpha, color=SCATTER_C, edgecolors="none",
                   rasterized=True, zorder=1)
        ax.plot([0, AXIS_MAX], [0, AXIS_MAX], "k--", lw=1.5, label="y = x", zorder=3)
        ax.plot(lx, ly, color=LOESS_C, lw=2.5, label="loess", zorder=4)
        ax.set_xlim(0, AXIS_MAX); ax.set_ylim(0, AXIS_MAX)
        ax.set_xlabel("True mutation age (generations)", fontsize=13)
        ax.set_ylabel("Midpoint estimate (generations)", fontsize=13)
        ax.legend(loc="upper left", fontsize=13, frameon=False)
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
        ax.spines[["top", "right"]].set_visible(False)
        if titled:
            ax.set_title("Bias in the midpoint age estimator (neutral coalescent simulation)")
            out = f"{out_prefix}.png"
        else:
            out = f"{out_prefix}.notitle.png"
        fig.tight_layout()
        fig.savefig(out, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"wrote {out}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rebuild", action="store_true",
                    help="recompute subsample + loess from the raw combined pickle")
    ap.add_argument("--out", default=os.path.join(HERE, "SuppFig_midpoint_bias_generations_100k"),
                    help="output path prefix (writes .png and .notitle.png)")
    ap.add_argument("--subsample", default=None,
                    help="comma-separated scatter subsample sizes drawn from the raw combined "
                         "pickle (e.g. 10000,50000); the loess is taken from the archived cache. "
                         "Writes <out>.sub<N>{,.notitle}.png per size.")
    a = ap.parse_args()

    if a.subsample:
        # larger scatter subsamples from the raw data, reusing the archived (all-mutation) loess
        _, _, lx, ly, _ = load_cache()
        print(f"loading {COMBINED} for larger subsamples ...")
        with open(COMBINED, "rb") as f:
            arr = np.asarray(pickle.load(f), dtype=np.int64)
        true_age, mid_age = arr[:, 0], arr[:, 1]
        rng = np.random.default_rng(SEED)
        for tok in a.subsample.split(","):
            k = int(tok)
            idx = rng.choice(len(true_age), size=k, replace=False)
            print(f"  subsample {k:,}")
            make_figure(true_age[idx], mid_age[idx], lx, ly, f"{a.out}.sub{k}")
        return

    ax_s, mx_s, lx, ly, xmax = rebuild_cache() if a.rebuild else load_cache()
    # Report the LAST downward crossing of y = x above 1 Ne (where the estimator turns
    # durably downward-biased). Near the origin the curve wiggles around y = x, so we ignore
    # crossings below 1 Ne and only report the meaningful high-age one.
    d = ly - lx
    cr = [i for i in np.where(np.diff(np.sign(d)) != 0)[0]
          if lx[i] > 10_000 and d[i] > 0 >= d[i + 1]]      # rising-x, +->- (downward) crossing
    if cr:
        i = cr[-1]
        xc = lx[i] + (lx[i+1] - lx[i]) * (0 - d[i]) / (d[i+1] - d[i])
        print(f"loess crosses y=x (downward) at {xc:,.0f} generations = {xc/10000:.2f} Ne")
    make_figure(ax_s, mx_s, lx, ly, a.out)


if __name__ == "__main__":
    main()
