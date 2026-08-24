#!/usr/bin/env python3
"""Methods accounting: how many candidate SNPs were actually tested for TSP status, and why not.

A candidate is counted as TESTED if its position appears in its region's `all_shared_metrics.txt`,
which means it survived chimpanzee recalling, was placed on the SINGER ARG, and had its coalescent
topology evaluated. Everything else is classified by the reason it was not:

  not_in_target_region                        outside every non-overlapping recall region
  not_tested__candidate_lost_in_chimp_recalling   region excluded (see check_candidate_presence.py)
  not_tested__region_not_analysed             SINGER metrics for the region do not exist
  not_tested__candidate_absent_from_ARG       region analysed, but this position is not in it

This is where the "tested N of 860 candidates" figure in the Methods comes from.

Inputs:  the SINGER working directory (candidate list, region list, exclusions, singer_trees/).
Output:  candidate_testing_summary.tsv, one row per candidate; breakdown printed.
"""
import os

import pandas as pd

# ---- inputs (edit) -------------------------------------------------------------------------
# The SINGER working directory: where ../../3_SINGER/ was run.
WORK = "../../3_SINGER"
OUT = "candidate_testing_summary.tsv"
# --------------------------------------------------------------------------------------------

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, WORK)

CAND = f"{WORK}/misc_files/merged_candidate_snps.chrom_pos.txt"
TARGET = f"{WORK}/misc_files/non_overlapping_region_250kb_ids.txt"
EXCL = f"{WORK}/misc_files/excluded_regions.txt"

candidates = set()
for line in open(CAND):
    p = line.split()
    if len(p) >= 2 and p[0] != "X":
        candidates.add((int(p[0]), int(p[1])))
candidates = sorted(candidates)

# target regions, as intervals per chromosome
regions = {}
for line in open(TARGET):
    r = line.strip()
    if not r:
        continue
    c, s, e = r.split("_")
    regions.setdefault(int(c), []).append((int(s), int(e), r))

excluded = {line.split()[0] for line in open(EXCL)
            if line.strip() and not line.startswith("#")}


def find_region(chrom, pos):
    for s, e, r in regions.get(chrom, []):
        if s <= pos <= e:
            return r
    return None


_tested_cache = {}


def tested_positions(region):
    """Positions in a region whose topology SINGER actually evaluated."""
    if region not in _tested_cache:
        f = f"{WORK}/singer_trees/{region}/{region}.all_shared_metrics.txt"
        positions = set()
        if os.path.exists(f):
            try:
                positions = set(pd.read_csv(f, sep="\t")["position"].astype(int))
            except Exception:
                pass
        _tested_cache[region] = positions
    return _tested_cache[region]


def has_metrics(region):
    return os.path.exists(
        f"{WORK}/singer_trees/{region}/{region}.all_shared_metrics.grouped.txt")


rows = []
for chrom, pos in candidates:
    region = find_region(chrom, pos)
    if region is None:
        category = "not_in_target_region"
    elif region in excluded:
        category = "not_tested__candidate_lost_in_chimp_recalling"
    elif not has_metrics(region):
        category = "not_tested__region_not_analysed"
    elif pos in tested_positions(region):
        category = "tested_for_TSP"
    else:
        category = "not_tested__candidate_absent_from_ARG"
    rows.append({"chrom": chrom, "pos": pos, "region": region, "category": category})

df = pd.DataFrame(rows)
df.to_csv(os.path.join(HERE, OUT), sep="\t", index=False)

print(f"TOTAL candidate SNPs: {len(df)}\n")
for k, v in df["category"].value_counts().items():
    print(f"  {v:>4}  {k}  ({100 * v / len(df):.1f}%)")
n_target = sum(len(v) for v in regions.values())
print(f"\nregions: {n_target} target, {len(excluded)} excluded")
print(f"wrote {OUT}")
