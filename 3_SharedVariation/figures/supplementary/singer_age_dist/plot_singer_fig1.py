#!/usr/bin/env python3
"""SINGER analogue of main Figure 1A: percentage of SNPs older than age,
shared vs all common 1000GP (MAF>=5%), SINGER posterior ages (100 ARG samples).
Styling matches Figure 1 Panel A."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib.patches import Patch

BASE = "../../../../3_SharedVariation/figures/supplementary/singer_age_dist"
plain_log_fmt = FuncFormatter(lambda y, _: f'{y:g}')

def gradient_vspan(ax, core, fade, color='orange', alpha_max=0.3, n=200, zorder=0):
    lo, hi = core; flo, fhi = lo - fade, hi + fade
    edges = np.linspace(flo, fhi, n + 1)
    for x0, x1 in zip(edges[:-1], edges[1:]):
        xm = 0.5 * (x0 + x1)
        if   xm < lo: a = alpha_max * (xm - flo) / (lo - flo)
        elif xm > hi: a = alpha_max * (fhi - xm) / (fhi - hi)
        else:         a = alpha_max
        ax.axvspan(x0, x1, color=color, alpha=max(a, 0.0), lw=0, zorder=zorder)
    return Patch(facecolor=color, alpha=alpha_max)

sh = np.loadtxt(f"{BASE}/outputs/singer_age_dist.shared.txt", skiprows=1)
al = np.loadtxt(f"{BASE}/outputs/singer_age_dist.all_common.txt", skiprows=1)
# cols: age_gen, age_mya, expected_snps, pct_older
sh_n = int(round(sh[0, 2])); al_n = int(round(al[0, 2]))

def make(ymin, suffix, xmax=10):
    fig, ax = plt.subplots(figsize=(5.2, 4.2), dpi=400)
    h_c_split = [5, 6]
    ax.plot(sh[:, 1], sh[:, 3], label=f"Shared SNPs (n={sh_n:,})", color='#1f77b4')
    ax.plot(al[:, 1], al[:, 3], label=f"1000GP (MAF ≥ 5%) (n={al_n:,})", color='red')
    ax.set_xlabel("Age (millions of years)", fontweight="bold")
    ax.set_ylabel("Percentage of SNPs (%)", fontweight="bold")
    ax.set_yscale('log')
    ax.set_xlim(0, xmax)
    ax.set_ylim(ymin, 100)
    ax.yaxis.set_major_formatter(plain_log_fmt)
    # integer x-ticks (e.g. 0,5,10,15,20 -> "5","10", not "5.0")
    step = 5 if xmax > 10 else 2
    ax.set_xticks(list(range(0, int(xmax) + 1, step)))
    ax.set_box_aspect(1)   # square plotting box
    hc_patch = gradient_vspan(ax, core=(h_c_split[0], h_c_split[1]), fade=0.5)
    _h, _l = ax.get_legend_handles_labels()
    _h.append(hc_patch); _l.append('Human-Chimp split')
    ax.legend(_h, _l, loc="lower left", frameon=False)
    plt.tight_layout()
    out = f"{BASE}/figures/SuppFig_singer_age_dist{suffix}.png"
    plt.savefig(out, dpi=400, bbox_inches='tight'); plt.close(fig)
    print("saved", out)

make(1e-2, "")                          # full range (matches Fig 1A axes)
make(1.0, ".ymin1")                     # y-axis truncated to start at 1%
make(1.0, ".ymin1_xmax20", xmax=20)     # y from 1%, x out to 20 Mya
