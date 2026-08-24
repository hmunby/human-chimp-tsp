#!/usr/bin/env python3
"""Explain each excluded region: why its candidate SNP did not survive chimpanzee recalling.

`check_candidate_presence.py` says which regions lost their candidate; this says why. For every
excluded region, the candidate position is looked up in the recall pipeline's per-chromosome
pass/fail position files and given one of:

  PASS(called+passed)              present and passing (so the loss happened later, e.g. in
                                   merging or phasing)
  FAILED_QC(LowDP/LowGeno/LowMQ)   called in chimpanzee but filtered out
  NEVER_CALLED_in_chimp            no chimpanzee call at that position at all
  no_verdict_file                  the recall workdir has no position files for that chromosome

Inputs:  the SINGER working directory (candidate list, exclusions) and the chimpanzee recall
         working directory from ../../2_ChimpRecall/, which holds the per-chromosome
         filtered/{chrom}/filtered.{chrom}_{pass,fail}.pos files.
Output:  lost_exclusions_summary.tsv, one row per excluded region per candidate.
"""
import os

import pandas as pd

# ---- inputs (edit) -------------------------------------------------------------------------
WORK = "../../3_SINGER"          # SINGER working directory
RECALL = "../../2_ChimpRecall"   # chimpanzee recall working directory

OUT = "lost_exclusions_summary.tsv"
# --------------------------------------------------------------------------------------------

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, WORK)
RECALL = os.path.join(HERE, RECALL)

# Every excluded region is a candidate-loss region: its candidate SNP either was never called in
# chimpanzee, failed QC, or is present only through a cross-chromosome position collision. A few
# also carry a stale SINGER tree, but they are genuine losses, so all are reported.
excluded = [line.split()[0] for line in open(f"{WORK}/misc_files/excluded_regions.txt")
            if line.strip() and not line.startswith("#")]

candidates = {}
for line in open(f"{WORK}/misc_files/merged_candidate_snps.chrom_pos.txt"):
    p = line.split()
    if len(p) >= 2 and p[0] != "X":
        candidates.setdefault(p[0], set()).add(int(p[1]))


def verdict(chrom, pos):
    def load(suffix):
        f = f"{RECALL}/filtered/{chrom}/filtered.{chrom}_{suffix}.pos"
        if not os.path.exists(f):
            return None
        positions = set()
        for x in open(f):
            x = x.strip()
            if not x:
                continue
            coord = x.split(":")[-1]        # chrN:pos  or  chrN:start-end
            for part in coord.split("-"):
                try:
                    positions.add(int(part))
                except ValueError:
                    pass
        return positions

    passed, failed = load("pass"), load("fail")
    if passed is None:
        return "no_verdict_file"
    if pos in passed:
        return "PASS(called+passed)"
    if failed and pos in failed:
        return "FAILED_QC(LowDP/LowGeno/LowMQ)"
    return "NEVER_CALLED_in_chimp"


rows = []
for region in sorted(excluded, key=lambda x: (int(x.split("_")[0]), int(x.split("_")[1]))):
    chrom, start, end = region.split("_")
    in_region = sorted(x for x in candidates.get(chrom, ()) if int(start) <= x <= int(end))
    for pos in (in_region or [None]):
        v = verdict(chrom, pos) if pos else "no_candidate_in_region"
        rows.append({"region": region, "chrom": int(chrom), "candidate_pos_hg19": pos,
                     "chimp_verdict": v})

df = pd.DataFrame(rows)
df.to_csv(os.path.join(HERE, OUT), sep="\t", index=False)
print(f"wrote {OUT}: {df['region'].nunique()} lost-candidate regions "
      f"({len(df)} candidate rows)")
print("\nverdict counts:")
print(df["chimp_verdict"].value_counts().to_string())
print("\nby chromosome:", dict(sorted(df.groupby("chrom")["region"].nunique().items())))
