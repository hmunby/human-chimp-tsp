#!/usr/bin/env python3
"""Genetic (recombination-map) span of each TSP region, for Table 1.

A region's physical span says little on its own: 6 kb in a recombination hotspot is a very
different object from 6 kb in a cold region. Converting each region's endpoints to centimorgans
gives the span in units of recombination, which is what determines whether the TSPs inside it can
plausibly have been carried as one haplotype since the human-chimp split.

Both endpoints are linearly interpolated onto the map, and the span is the difference. Two maps
are used so the result does not rest on one:

    HapMap II (GRCh37)          the primary map, and the source of the mean rate reported
    1000 Genomes OMNI           an independent interpolated map, as a sensitivity check

Regions with a single TSP have zero span by construction; the multi-TSP regions are the
informative ones.

Coordinates are hg19/b37, matching both maps and the clustering.

Inputs:  tsp_clusters_10kb.pooled.tsv from aggregate_seeds.py, and (for the gene label only)
         ../3_GeneAssignment/region_gene_assignments.pooled.tsv; the two genetic maps.
Outputs: tsp_region_genetic_span.tsv             every region, both maps
         tsp_region_genetic_span_comparison.tsv  the same, formatted for the manuscript
"""
import argparse
import os

import numpy as np
import pandas as pd

# ---- inputs (edit) -------------------------------------------------------------------------
CLUSTERS = "tsp_clusters_10kb.pooled.tsv"
GENE_ASSIGNMENTS = "../3_GeneAssignment/region_gene_assignments.pooled.tsv"

# Genetic maps, GRCh37. Not bundled; from https://github.com/joepickrell/1000-genomes-genetic-maps
HAPMAP_TMPL = "resources/1000-genomes-genetic-maps/hapmapII/genetic_map_GRCh37_chr{chrom}.txt"
OMNI_TMPL = ("resources/1000-genomes-genetic-maps/interpolated_OMNI/"
             "chr{chrom}.OMNI.interpolated_genetic_map.gz")

OUT = "tsp_region_genetic_span.tsv"
OUT_COMPARISON = "tsp_region_genetic_span_comparison.tsv"

# Regions that within_region_ld.py showed are NOT a single linked block, and so are reported as
# their separate LD blocks in the manuscript table. Reporting one span across a region that holds
# two unlinked blocks would describe a distance no observed haplotype actually spans.
#
# IGFBP7 is the only such region. In the reported (pooled) set its six TSPs form two blocks --
# two at 57,918,296-57,918,492 and four at 57,919,221-57,919,705 -- with minimum pairwise
# r^2 of 0.024 (AFR), 0.003 (EUR) and 0.038 (YRI) across the split, against a mean of ~0.48.
# The block boundaries are unchanged from the draft set; only the second block's TSP count is.
#
#   region id -> [(label, start, end, n_tsps), ...]
LD_SUBGROUPS = {
    "4_57918296_57919705": [
        ("IGFBP7-A", 57918296, 57918492, 2),
        ("IGFBP7-B", 57919221, 57919705, 4),
    ],
}
# --------------------------------------------------------------------------------------------

HERE = os.path.dirname(os.path.abspath(__file__))


def path(p):
    return p if os.path.isabs(p) else os.path.join(HERE, p)


class GeneticMap:
    """Linear interpolation of position -> cM, loading each chromosome on first use."""

    def __init__(self, template, reader):
        self.template = template
        self.reader = reader
        self._cache = {}

    def _chrom(self, chrom):
        if chrom not in self._cache:
            self._cache[chrom] = self.reader(path(self.template.format(chrom=chrom)))
        return self._cache[chrom]

    def cM(self, chrom, pos):
        positions, cm = self._chrom(chrom)
        return float(np.interp(pos, positions, cm))

    def span(self, chrom, start, end):
        return self.cM(chrom, end) - self.cM(chrom, start)


def read_hapmap(f):
    m = pd.read_csv(f, sep="\t")
    return m["Position(bp)"].to_numpy(), m["Map(cM)"].to_numpy()


def read_omni(f):
    m = pd.read_csv(f, sep=r"\s+", header=None, names=["rsid", "pos", "cM"], compression="gzip")
    return m["pos"].to_numpy(), m["cM"].to_numpy()


# The values above are the published defaults; each can be overridden on the command line so
# the same measurement can be taken on an alternative TSP set (the Ne-sensitivity comparison
# uses this). --no-ld-subgroups turns off the IGFBP7 split, which is a fact about the published
# set and does not carry over to a set whose regions differ.
_ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
_ap.add_argument("--clusters", default=CLUSTERS)
_ap.add_argument("--genes", default=GENE_ASSIGNMENTS)
_ap.add_argument("--out", default=OUT)
_ap.add_argument("--out-comparison", default=OUT_COMPARISON)
_ap.add_argument("--hapmap", default=HAPMAP_TMPL, help="GRCh37 HapMap II map, with {chrom}")
_ap.add_argument("--omni", default=OMNI_TMPL, help="GRCh37 OMNI map, with {chrom}")
_ap.add_argument("--no-ld-subgroups", action="store_true",
                 help="do not split any region into LD sub-blocks")
_a = _ap.parse_args()
CLUSTERS, GENE_ASSIGNMENTS = _a.clusters, _a.genes
OUT, OUT_COMPARISON = _a.out, _a.out_comparison
HAPMAP_TMPL, OMNI_TMPL = _a.hapmap, _a.omni
if _a.no_ld_subgroups:
    LD_SUBGROUPS = {}

clusters = pd.read_csv(path(CLUSTERS), sep="\t")
clusters["region"] = (clusters["chrom"].astype(str) + "_"
                      + clusters["start"].astype(str) + "_" + clusters["end"].astype(str))

genes = {}
if os.path.exists(path(GENE_ASSIGNMENTS)):
    ga = pd.read_csv(path(GENE_ASSIGNMENTS), sep="\t")
    genes = dict(zip(ga["region"], ga["assigned_gene"]))
else:
    print(f"note: {GENE_ASSIGNMENTS} not found, gene column left blank")

hapmap = GeneticMap(HAPMAP_TMPL, read_hapmap)
omni = GeneticMap(OMNI_TMPL, read_omni)

rows = []
for c in clusters.itertuples():
    chrom = str(c.chrom)
    span_bp = c.end - c.start
    span_cM = hapmap.span(chrom, c.start, c.end)
    rows.append({
        "region": c.region,
        "gene": genes.get(c.region, ""),
        "n_tsps": c.n_tsps,
        "span_bp": span_bp,
        "span_cM": span_cM,
        "mean_rate_cM_Mb": (span_cM / (span_bp / 1e6)) if span_bp > 0 else np.nan,
        "span_cM_OMNI": omni.span(chrom, c.start, c.end),
    })

span = pd.DataFrame(rows)
span.to_csv(path(OUT), sep="\t", index=False)


def comparison_row(label, chrom, start, end, n_tsps):
    """One manuscript-table row: spans in 10^-3 cM, region as chr:start-end."""
    span_bp = end - start
    hm = hapmap.span(chrom, start, end) * 1e3
    om = omni.span(chrom, start, end) * 1e3
    return {
        "Region (hg19)": f"chr{chrom}:{start:,}-{end:,}",
        "Gene": label,
        "TSPs": n_tsps,
        "Span (bp)": span_bp,
        "Span_1e-3cM_HapMapII": round(hm, 3),
        "Span_1e-3cM_OMNI": round(om, 3),
        "MeanRate_cM_Mb_HapMapII": round((hm / 1e3) / (span_bp / 1e6), 2) if span_bp > 0 else np.nan,
    }


rows = []
for c in clusters.itertuples():
    if c.region in LD_SUBGROUPS:
        for label, start, end, n in LD_SUBGROUPS[c.region]:
            rows.append(comparison_row(label, str(c.chrom), start, end, n))
    else:
        rows.append(comparison_row(genes.get(c.region, ""), str(c.chrom),
                                   c.start, c.end, c.n_tsps))

comparison = pd.DataFrame(rows).sort_values(["TSPs", "Span (bp)"], ascending=False)
comparison.to_csv(path(OUT_COMPARISON), sep="\t", index=False)

multi = comparison[comparison["TSPs"] > 1]
print("Genetic span of multi-TSP regions (10^-3 cM), HapMap II vs 1000G OMNI, GRCh37:\n")
print(f"  {'Gene':>9} {'TSPs':>4} {'Span (bp)':>9} {'HapMapII':>9} {'OMNI':>8} {'cM/Mb':>7}")
for _, r in multi.iterrows():
    print(f"  {str(r['Gene']):>9} {r['TSPs']:>4} {r['Span (bp)']:>9,}"
          f" {r['Span_1e-3cM_HapMapII']:>9.3f} {r['Span_1e-3cM_OMNI']:>8.3f}"
          f" {r['MeanRate_cM_Mb_HapMapII']:>7.2f}")

n_split = sum(len(v) for v in LD_SUBGROUPS.values()) - len(LD_SUBGROUPS)
print(f"\n  total over the {len(multi)} multi-TSP rows "
      f"({len(multi) - n_split} regions, IGFBP7 split into its two LD blocks): "
      f"HapMapII={multi['Span_1e-3cM_HapMapII'].sum():.2f}, "
      f"OMNI={multi['Span_1e-3cM_OMNI'].sum():.2f} (x10^-3 cM)")
print(f"  largest span under either map: "
      f"{max(multi['Span_1e-3cM_HapMapII'].max(), multi['Span_1e-3cM_OMNI'].max())/1e3:.5f} cM")
print(f"  largest disagreement between maps: "
      f"{(multi['Span_1e-3cM_HapMapII'] - multi['Span_1e-3cM_OMNI']).abs().max():.2f} x10^-3 cM")
print(f"\nwrote {OUT} ({len(span)} regions) and {OUT_COMPARISON} ({len(comparison)} rows)")
