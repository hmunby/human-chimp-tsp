#!/usr/bin/env python3
"""phyloP eCDF across SNP sets (Panel B of the combined variant-effects/phyloP figure).

Five tracks, drawn as empirical cumulative distributions of per-SNP phyloP score:

    Chimpanzee (AC >= 3)      black
    1000GP (MAF >= 5%)        grey     )
    All shared                blue     ) nested: 1000GP > all shared > old shared > TSP
    Old shared (>= 4 Mya)     purple   )
    TSP                       red      )

Run twice, once per phyloP track:
    phylop/       36-eutherian-mammal EPO, human sequence EXCLUDED  -> the main paper figure
    phylop_46way/ UCSC 46-way placental, human INCLUDED             -> comparison

Old shared and TSP are subsets of the shared SNPs, so their scores are looked up in the
shared table rather than scored separately.

Inputs:  phylop*/ per-chromosome score tables from ./Snakefile; the two SNP-set position
         lists in misc_files/.
Outputs: phylop_ecdf_{excl,incl}_humans.png (+ .notitle.png).

Importable: combined_variant_effects_phylop.py reuses _compute_tracks() and draw_ecdf().
"""
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))   # cwd-robust so this stays importable

# ---- inputs (edit) -------------------------------------------------------------------------
OLD_SHARED_FILE = os.path.join(HERE, "misc_files/old_shared_snps_860.hg19.chrom_pos.txt")
TSP_SNP_FILE = os.path.join(HERE, "misc_files/tsp_snps_pooled.hg19.txt")
# --------------------------------------------------------------------------------------------

# Shared SNP-set palette (nested gradient; TSP = bright red). Chimpanzee has no counterpart
# in Panel A of the combined figure, so it keeps black.
COLORS = {
    "chimp": "#000000",
    "tgp": "#7F7F7F",
    "allshared": "#1F77B4",
    "old4": "#9467BD",
    "tsp": "#FF0000",
}


def _load_chrom_pos(path):
    """Read a CHROM/POS table (TSP / old-shared) whose header contains CHROM and POS."""
    df = pd.read_csv(path, sep="\t", low_memory=False)
    df = df.rename(columns={c: c.upper() for c in df.columns})
    df["CHROM"] = pd.to_numeric(df["CHROM"], errors="coerce")
    df["POS"] = pd.to_numeric(df["POS"], errors="coerce")
    return df.dropna(subset=["CHROM", "POS"]).astype({"CHROM": int, "POS": int})


def _compute_tracks(phylop_dir):
    """Return {track name: phyloP array, NaNs dropped} for one phyloP source directory."""
    phylop_dir = os.path.join(HERE, phylop_dir)
    shared_ph = {c: pd.read_csv(f"{phylop_dir}/shared_snps_phylop_chr{c}.txt", sep="\t")
                 for c in range(1, 23)}

    def flat(prefix):
        out = []
        for c in range(1, 23):
            out += pd.read_csv(f"{phylop_dir}/{prefix}_chr{c}.txt", sep="\t")["PHYLOP_SCORE"].tolist()
        return np.asarray(out, dtype=float)

    def subset_of_shared(chrom_pos):
        out = []
        for c in range(1, 23):
            sub = chrom_pos[chrom_pos["CHROM"] == c]
            if len(sub):
                m = pd.merge(sub[["CHROM", "POS"]], shared_ph[c], on=["CHROM", "POS"], how="inner")
                out += m["PHYLOP_SCORE"].tolist()
        return np.asarray(out, dtype=float)

    tracks = {
        "chimp": flat("chimp_snps_phylop"),
        "tgp": flat("tgp_maf5_snps_phylop"),
        "allshared": np.concatenate(
            [shared_ph[c]["PHYLOP_SCORE"].to_numpy(dtype=float) for c in range(1, 23)]),
        "old4": subset_of_shared(_load_chrom_pos(OLD_SHARED_FILE)),
        "tsp": subset_of_shared(_load_chrom_pos(TSP_SNP_FILE)),
    }
    return {k: v[~np.isnan(v)] for k, v in tracks.items()}


def draw_ecdf(ax, tracks, legend_fs="small", label_fs=None, tick_fs=None,
              loc="upper left", handlelength=1.2):
    """Draw the five-track eCDF onto `ax`. No title; every legend entry carries its n.

    TSP is drawn last (and thicker) so the red line stays on top, but the legend order is
    set explicitly so it does not follow draw order: Chimpanzee first, then the same set
    order as Panel A of the combined figure. Font sizes default to the standalone values;
    the combined-figure driver passes larger ones.
    """
    n = {k: len(v) for k, v in tracks.items()}
    ax.axvline(0, color="grey", ls="--", lw=1, alpha=0.5, zorder=0)   # neutral reference

    # ax.ecdf returns a single Line2D (not a container), so do not unpack.
    h_chimp = ax.ecdf(tracks["chimp"], color=COLORS["chimp"])
    h_tgp = ax.ecdf(tracks["tgp"], color=COLORS["tgp"])
    h_all = ax.ecdf(tracks["allshared"], color=COLORS["allshared"])
    h_old = ax.ecdf(tracks["old4"], color=COLORS["old4"])
    h_tsp = ax.ecdf(tracks["tsp"], color=COLORS["tsp"], linewidth=2)

    labels = [
        f"Chimpanzee (AC ≥ 3) (n={n['chimp']:,})",
        f"1000GP (MAF ≥ 5%) (n={n['tgp']:,})",
        f"All shared (n={n['allshared']:,})",
        f"Old shared (≥4 Mya) (n={n['old4']:,})",
        f"TSP (n={n['tsp']:,})",
    ]
    ax.set_xlabel("Phylop Score", fontweight="bold", fontsize=label_fs)
    ax.set_ylabel("Cumulative Density", fontweight="bold", fontsize=label_fs)
    ax.set_xlim(-5, 4)
    if tick_fs is not None:
        ax.tick_params(labelsize=tick_fs)
    ax.legend([h_chimp, h_tgp, h_all, h_old, h_tsp], labels, loc=loc, frameon=False,
              fontsize=legend_fs, handlelength=handlelength, handletextpad=0.4,
              labelspacing=0.3, borderaxespad=0.3)
    return ax


OUT_DIR = None   # set by --out-dir; None means alongside this script


def make_plot(phylop_dir, version, out_png):
    tracks = _compute_tracks(phylop_dir)
    print(f"\n=== {version} ({phylop_dir}) ===")
    for k, v in tracks.items():
        print(f"  {k:10s} n={len(v):>9d}  mean={np.mean(v):+.3f}")

    fig, ax = plt.subplots(figsize=(5, 4), dpi=400)
    draw_ecdf(ax, tracks)
    fig.tight_layout()
    out_png = os.path.join(OUT_DIR or HERE, out_png)
    fig.savefig(out_png.replace(".png", ".notitle.png"), bbox_inches="tight")
    # the titled version keeps incl/excl humans in the title to tell the two files apart
    ax.set_title(f"Cumulative Density of Phylop Scores ({version})")
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {out_png} (+ .notitle)")


if __name__ == "__main__":
    import argparse
    _ap = argparse.ArgumentParser()
    _ap.add_argument("--tsp-file", help="alternative TSP set; defaults to the reported pooled set")
    _ap.add_argument("--out-dir", help="write the figures here instead of alongside the script")
    _a, _rest = _ap.parse_known_args()
    if _a.tsp_file:
        TSP_SNP_FILE = _a.tsp_file
    if _a.out_dir:
        import os as _os
        _os.makedirs(_a.out_dir, exist_ok=True)
        OUT_DIR = _a.out_dir      # output only; HERE still locates phylop/ and misc_files/

    make_plot("phylop", "excl. humans", "phylop_ecdf_excl_humans.png")
    make_plot("phylop_46way", "incl. humans", "phylop_ecdf_incl_humans.png")
