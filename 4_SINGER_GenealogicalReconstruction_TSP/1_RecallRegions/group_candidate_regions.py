#!/usr/bin/env python3
"""Group candidate TSP SNPs into recall windows for stage 4.

Each candidate SNP (from ../../3_SharedVariation/2_CandidateTSPs/) is expanded to a +/- 250 kb
window; overlapping windows on a chromosome are merged into non-overlapping recall regions. These
regions decide where chimp variants are recalled (../2_ChimpRecall/) and where SINGER reconstructs
genealogies (../3_SINGER/).

NOTE: this grouping is only to choose recall windows. It is NOT the region grouping of the final
TSPs we report -- those are derived later from the SINGER output.

Outputs (consumed by ../2_ChimpRecall/):
  regions_chr{chrom}.list        GATK interval list, one candidate SNP per line (padded +/-250kb by
                                 HaplotypeCaller's -ip 250000)
  regions_chr{chrom}_250kb.txt   merged +/- 250 kb windows (CHROM<TAB>START<TAB>END)

The exact region files used in the paper are bundled with the recalling pipeline
(../2_ChimpRecall/misc_files/regions_chr*.list, ../2_ChimpRecall/region_files/regions_chr*_250kb*.txt);
use those for exact reproduction. This script documents the grouping logic.
"""
from collections import defaultdict

CANDIDATES = "../../3_SharedVariation/2_CandidateTSPs/candidate_tsps.chrom_pos.txt"
PAD = 250_000

by_chrom = defaultdict(list)
with open(CANDIDATES) as f:
    for line in f:
        c, p = line.split()[:2]
        by_chrom[c].append(int(p))

for c, positions in by_chrom.items():
    positions = sorted(positions)
    # GATK interval list: one candidate SNP per line (HaplotypeCaller pads +/-250kb via -ip)
    with open(f"regions_chr{c}.list", "w") as f:
        for p in positions:
            f.write(f"{c}:{p}\n")
    # merge +/- 250 kb windows on this chromosome
    merged = []
    for p in positions:
        s, e = max(1, p - PAD), p + PAD
        if merged and s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    with open(f"regions_chr{c}_250kb.txt", "w") as f:
        for s, e in merged:
            f.write(f"{c}\t{s}\t{e}\n")
