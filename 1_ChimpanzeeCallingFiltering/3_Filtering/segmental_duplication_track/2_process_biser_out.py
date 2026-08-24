#!/usr/bin/env python3
"""Segmental-duplication track, step 2: extract panTro6 SD regions that are resolved in T2T.

From the BISER output (step 1), keep only the cross-assembly pairs "T2T:panTro6" -- segmental
duplications that BISER anchored in panTro6 and that are resolved (single-copy) in the T2T mPanTro3
assembly. Emit their panTro6 coordinates as BED, with and without the "chr" prefix. Overlapping
intervals are merged in step 3.
"""
import os
import pandas as pd

BISER_OUT = "out/t2t_pantro6_seg_dups.txt"
OUT_DIR = "processed_bed"

# assembly labels exactly as they appear in BISER's "assemblies" column
T2T = "mPanTro3.hap1.cur.20231122.masked"
PANTRO6 = "panTro6"

# panTro6 chromosome names (chimp chr2 is split into 2A / 2B)
CHROMS = ["chr1", "chr2A", "chr2B"] + [f"chr{i}" for i in range(3, 23)] + ["chrX", "chrY"]

COLS = ["chrom_1", "start_1", "end_1", "chrom_2", "start_2", "end_2", "assemblies", "score",
        "strand_1", "strand_2", "len_1", "len_2", "cigar", "extra"]

biser = pd.read_csv(BISER_OUT, sep="\t", header=None, names=COLS)

# cross-assembly pairs anchored on the panTro6 side (chrom_2 = panTro6 coordinates)
resolved = biser[biser["assemblies"] == f"{T2T}:{PANTRO6}"]
resolved = resolved[resolved["chrom_2"].isin(CHROMS)]          # drop unplaced contigs
resolved = resolved[["chrom_2", "start_2", "end_2"]]

# order by chromosome then position, dedupe
resolved["chrom_2"] = pd.Categorical(resolved["chrom_2"], categories=CHROMS, ordered=True)
resolved = resolved.sort_values(["chrom_2", "start_2"]).drop_duplicates()

os.makedirs(OUT_DIR, exist_ok=True)
resolved.to_csv(f"{OUT_DIR}/t2t_pantro6_seg_dups_resolved.bed",
                sep="\t", header=False, index=False)

no_chr = resolved.copy()
no_chr["chrom_2"] = no_chr["chrom_2"].astype(str).str.replace("chr", "", regex=False)
no_chr.to_csv(f"{OUT_DIR}/t2t_pantro6_seg_dups_resolved_no_chr.bed",
              sep="\t", header=False, index=False)
