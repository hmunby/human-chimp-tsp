#!/usr/bin/env python3
"""Combined two-panel figure: variant effects (A) and phyloP eCDF (B).

Reuses the computations and plotters from variant_effects.py and phylop_ecdf.py so the panels
stay identical to the standalone figures and share the SNP-set palette. One combined figure
per phyloP source (humans excluded / included), in side-by-side and stacked layouts.

Outputs: combined_variant_effects_phylop.{excl,incl}_humans{,.vertical}.png
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import phylop_ecdf as ph          # exposes _compute_tracks + draw_ecdf
import variant_effects as ve      # builds `counts` at import (no plotting)

# Panel A data (merged categories)
pctA, seA, nsA = ve.effect_table(ve.counts, ve.ORDER_MERGED, merge=ve.MERGE_GROUPS)

# Enlarged fonts for the combined figure; the standalone scripts keep their own defaults.
TICK_FS, LEGEND_FS, LABEL_FS, PANEL_FS, MARKER_MS = 14, 13, 16, 19, 6
LEGEND_FS_B = 10   # Panel B legend smaller, with shorter handles, so it clears the eCDF curves

LAYOUTS = {
    "": dict(shape=(1, 2), figsize=(14, 7.2), labxy=(-0.18, 1.02)),           # side by side
    ".vertical": dict(shape=(2, 1), figsize=(7.5, 15), labxy=(-0.22, 1.02)),  # stacked
}

for phylop_dir, tag in [("phylop", "excl_humans"), ("phylop_46way", "incl_humans")]:
    tracks = ph._compute_tracks(phylop_dir)   # compute once, reuse across layouts
    for suffix, layout in LAYOUTS.items():
        fig, (axA, axB) = plt.subplots(*layout["shape"], figsize=layout["figsize"], dpi=300)

        ve.plot_effects(pctA, seA, nsA, ve.ORDER_MERGED, title=None, ax=axA,
                        tick_fs=TICK_FS, legend_fs=LEGEND_FS, label_fs=LABEL_FS,
                        marker_ms=MARKER_MS)
        ph.draw_ecdf(axB, tracks, legend_fs=LEGEND_FS_B, label_fs=LABEL_FS,
                     tick_fs=TICK_FS, handlelength=0.7)

        axA.set_box_aspect(1)
        axB.set_box_aspect(1)

        lx, ly = layout["labxy"]
        for ax, lab in [(axA, "A"), (axB, "B")]:
            ax.text(lx, ly, lab, transform=ax.transAxes, fontsize=PANEL_FS,
                    fontweight="bold", va="bottom", ha="right")

        fig.tight_layout()
        out = os.path.join(HERE, f"combined_variant_effects_phylop.{tag}{suffix}.png")
        fig.savefig(out, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"saved {out}")
