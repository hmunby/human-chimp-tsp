#!/usr/bin/env python3
"""Variant-effect (snpEff consequence) distribution across SNP sets.

Panel A of the combined variant-effects/phyloP figure. Four nested sets:
1000GP (MAF >= 5%) > all shared > old shared (>= 4 Mya) > TSP.

MULTIPLE ANNOTATIONS PER SNP -> GREATEST CONSEQUENCE. snpEff reports one effect per
overlapping transcript. For each SNP we keep the effect ranking highest in CONSEQUENCE_ORDER
below, NOT snpEff's four-level IMPACT. The category is therefore re-derived here from the raw
multi-annotation snpEff file rather than taken from the upstream IMPACT-collapsed summaries
(this changes ~4% of shared-SNP categories).

Categories are ordered least -> greatest consequence. Coding and regulatory categories are
ordered by the constraint estimate `est` from Agarwal et al. Fig 2B (synonymous ~ neutral 1.0
down to LOF ~ 0.73); the non-coding categories Agarwal does not cover (intergenic, intron,
intragenic, non-coding exon) are placed at the low-consequence end.

Two figures are produced:
    .full   every effect kept separate
    .merged 5' UTR + 3' UTR -> UTR; upstream + downstream -> flanking

Old shared and TSP are subsets of the shared SNPs, so their categories are looked up in the
shared map, guaranteeing identical treatment. 1000GP uses a precomputed count cache built with
the same greatest-consequence rule over its (2.8 GB) snpEff file.

Importable: combined_variant_effects_phylop.py reuses `counts`, effect_table() and
plot_effects().
"""
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))   # cwd-robust so this stays importable

# ---- inputs (edit) -------------------------------------------------------------------------
# Raw snpEff annotations of the shared variants: CHROM, POS, ANN, LOF (no header). Produced by
# `bcftools query -f '%CHROM\t%POS\t%ANN\t%LOF\n'` on the stage-3 annotated shared VCF.
# ~27 MB, not bundled -- see README.
SHARED_SNPEFF = os.path.join(HERE, "resources/shared_variants.SD_CNV_DUP_filtered.snpeff.txt")

TGP_CACHE = os.path.join(HERE, "misc_files/tgp_annotation_simple_counts.txt")
OLD_FILE = os.path.join(HERE, "misc_files/old_shared_snps_860.hg19.chrom_pos.txt")
TSP_FILE = os.path.join(HERE, "misc_files/tsp_snps_pooled.hg19.txt")
AGARWAL = os.path.join(HERE, "misc_files/agarwal_et_al_cpg_fig2b.txt")
# --------------------------------------------------------------------------------------------

# Shared SNP-set palette (nested gradient; TSP = bright red).
COLORS = {"1000GP (MAF ≥ 5%)": "#7F7F7F", "All shared": "#1F77B4",
          "Old shared (≥4 Mya)": "#9467BD", "TSP": "#FF0000"}

# Consequence ordering (least -> most).
agarwal_order = pd.read_csv(AGARWAL, sep="\t").sort_values("est", ascending=False)["csq"].tolist()
LOW_END = ["intergenic", "intron", "intragenic", "non_coding_transcript_exon"]  # not in Agarwal
CONSEQUENCE_ORDER = LOW_END + agarwal_order
RANK = {c: i for i, c in enumerate(CONSEQUENCE_ORDER)}   # higher rank = greater consequence

# Group merges used only by the .merged figure
MERGE_GROUPS = {"5_prime_UTR": "UTR", "3_prime_UTR": "UTR",
                "upstream_gene": "flanking", "downstream_gene": "flanking"}
ORDER_MERGED = list(dict.fromkeys(MERGE_GROUPS.get(c, c) for c in CONSEQUENCE_ORDER))

LABEL_OVERRIDES = {"non_coding_transcript_exon": "non coding transcript"}


def cat_of_annotation(a):
    """One snpEff effect string -> simplified category."""
    if "stop_gained" in a or "start_lost" in a or "stop_lost" in a:
        return "LOF"
    if "missense_variant" in a:
        return "missense"
    if "stop_retained_variant" in a:
        return "synonymous"
    if "synonymous_variant" in a:
        return "synonymous"
    if "initiator_codon_variant" in a:
        return "missense"
    if "splice" in a:
        return "splice_region"
    if "upstream_gene_variant" in a:
        return "upstream_gene"
    if "downstream_gene_variant" in a:
        return "downstream_gene"
    if "5_prime_UTR" in a:
        return "5_prime_UTR"
    if "3_prime_UTR" in a:
        return "3_prime_UTR"
    if "non_coding_transcript_exon_variant" in a:
        return "non_coding_transcript_exon"
    if "intron_variant" in a:
        return "intron"
    if "intergenic_region" in a:
        return "intergenic"
    if "intragenic" in a:
        return "intragenic"
    return "other"


def best_category(ann_field, lof):
    """Among a SNP's annotations, the category with the greatest consequence."""
    if lof != ".":
        return "LOF"
    return max((cat_of_annotation(a) for a in str(ann_field).split(",")),
               key=lambda c: RANK.get(c, -1))


def effect_table(counts, order, merge=None):
    """Per-set percentages and binomial standard errors, optionally collapsing categories."""
    pct, se, ns = {}, {}, {}
    for name, s in counts.items():
        if merge:
            s = s.groupby(lambda c: merge.get(c, c)).sum()
        n = int(s.sum())
        ns[name] = n
        p = s.reindex(order, fill_value=0) / n
        pct[name] = p * 100
        se[name] = np.sqrt(p * (1 - p) / n) * 100
    return pct, se, ns


def plot_effects(pct, se, ns, order, title, savename=None, ax=None,
                 tick_fs=9, legend_fs=8, label_fs=None, marker_ms=4):
    """Draw the ordered dot plot (percent +/- SE, log x, greatest consequence at the top).

    With `ax`, draw onto it and return without saving (for embedding as a panel). Otherwise
    create a figure and save titled + .notitle versions to `savename`. Font sizes default to
    the standalone values; the combined-figure driver passes larger ones.
    """
    own = ax is None
    if own:
        fig, ax = plt.subplots(figsize=(7, 6), dpi=400)
    order_top_down = order[::-1]
    y = np.arange(len(order_top_down))
    for name, off in zip(pct, np.linspace(-0.24, 0.24, len(pct))):
        ax.errorbar(pct[name].reindex(order_top_down).values, y + off,
                    xerr=se[name].reindex(order_top_down).values,
                    fmt="o", ms=marker_ms, capsize=3, color=COLORS[name],
                    label=f"{name} (n={ns[name]:,})")
    ax.set_yticks(y)
    ax.set_yticklabels([LABEL_OVERRIDES.get(c, c.replace("_", " ")) for c in order_top_down],
                       fontsize=tick_fs, rotation=45, ha="right", va="center",
                       rotation_mode="anchor")
    ax.set_xscale("log")
    ax.tick_params(axis="x", labelsize=tick_fs)
    ax.set_xlabel("Percentage of SNPs (%)", fontweight="bold", fontsize=label_fs)
    ax.set_ylabel("Variant effect", fontweight="bold", fontsize=label_fs)
    ax.invert_yaxis()
    ax.legend(frameon=False, fontsize=legend_fs, loc="lower left")
    ax.grid(axis="x", which="both", ls=":", alpha=0.3)
    if not own:
        return ax
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, savename.replace(".png", ".notitle.png")), dpi=300)
    ax.set_title(title, fontweight="bold")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, savename), dpi=300)
    plt.close(fig)
    return ax


_SHARED_MAP = None


def shared_category_map():
    """(CHROM, POS) -> greatest-consequence category, for every shared SNP.

    Cached: the snpEff table is 27 MB and re-deriving it per TSP set dominates the cost
    when several sets are scored in one process.
    """
    global _SHARED_MAP
    if _SHARED_MAP is None:
        raw = pd.read_csv(SHARED_SNPEFF, sep="\t", header=None,
                          names=["CHROM", "POS", "ANN", "LOF"])
        raw["cat"] = [best_category(a, l) for a, l in zip(raw["ANN"], raw["LOF"])]
        _SHARED_MAP = raw.set_index(["CHROM", "POS"])["cat"]
    return _SHARED_MAP


def subset_counts(path, shared_map=None):
    """Category counts for a subset of the shared SNPs, given as a CHROM/POS table."""
    shared_map = shared_category_map() if shared_map is None else shared_map
    d = pd.read_csv(path, sep="\t", low_memory=False)
    cats = shared_map.reindex(list(zip(d["CHROM"], d["POS"])))
    missing = cats.isna().sum()
    if missing:
        print(f"  WARNING: {missing} of {len(d)} positions in {os.path.basename(path)} "
              "not found in the shared snpEff table")
    return cats.dropna().value_counts()


def _build_counts(tsp_file=None, old_file=None):
    """Per-set category counts. Runs at import so the combined figure can reuse `counts`.

    The two subset paths are arguments so an alternative TSP set can be scored through
    exactly this code (the Ne-sensitivity comparison does that).
    """
    shared_map = shared_category_map()
    return {
        "1000GP (MAF ≥ 5%)": pd.read_csv(TGP_CACHE, sep="\t").set_index("Annotation_Simple")["count"],
        "All shared": shared_map.value_counts(),
        "Old shared (≥4 Mya)": subset_counts(old_file or OLD_FILE, shared_map),
        "TSP": subset_counts(tsp_file or TSP_FILE, shared_map),
    }


counts = _build_counts()

_unknown = set().union(*[set(s.index) for s in counts.values()]) - set(CONSEQUENCE_ORDER)
if _unknown:
    print("WARNING: categories not in CONSEQUENCE_ORDER:", _unknown)
print("Set sizes (N):", {k: int(v.sum()) for k, v in counts.items()})


if __name__ == "__main__":
    import argparse
    _ap = argparse.ArgumentParser(description=__doc__.split("\n")[0] if __doc__ else None)
    _ap.add_argument("--tsp-file", help="alternative TSP set (CHROM/POS table); "
                                        "defaults to the reported pooled set")
    _ap.add_argument("--old-file", help="alternative old-shared set")
    _a = _ap.parse_args()
    if _a.tsp_file or _a.old_file:
        counts = _build_counts(_a.tsp_file, _a.old_file)
        print("Set sizes (N):", {k: int(v.sum()) for k, v in counts.items()})

    pct, se, ns = effect_table(counts, CONSEQUENCE_ORDER)
    plot_effects(pct, se, ns, CONSEQUENCE_ORDER,
                 "Variant effects — full categories", "variant_effects_by_set.full.png")

    pct, se, ns = effect_table(counts, ORDER_MERGED, merge=MERGE_GROUPS)
    plot_effects(pct, se, ns, ORDER_MERGED,
                 "Variant effects", "variant_effects_by_set.merged.png")
