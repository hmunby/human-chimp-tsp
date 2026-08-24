#!/usr/bin/env python3
"""Drop lifted segments whose length changed materially in the panTro6 -> hg38 liftOver.

Each minimal-decode segment is 1 kb in panTro6. A lift that returns something far from 1 kb
means the segment straddled a rearrangement or an indel-rich region, so its posterior no longer
describes a single comparable interval. Only 950-1050 bp segments are kept.

Writes the combined bed plus one bed per autosome (the form ../../2_UncoalescedSegments/ reads).

Usage: clean_lifted_minimal_decode.py <combined_hg38.bed> <out_dir> <popsam>
"""
import os
import sys

import pandas as pd

MIN_LEN, MAX_LEN = 950, 1050


def main(infile, out_dir, popsam):
    os.makedirs(out_dir, exist_ok=True)

    df = pd.read_csv(infile, sep="\t", header=None,
                     names=["chrom", "start", "end", "posterior"])
    length = df["end"] - df["start"]
    df = df[(length > MIN_LEN) & (length < MAX_LEN)]

    df.to_csv(f"{out_dir}/{popsam}_minimal_decode_hg38.bed",
              sep="\t", header=False, index=False)
    for chrom in range(1, 23):
        df[df["chrom"] == f"chr{chrom}"].to_csv(
            f"{out_dir}/{popsam}_chr{chrom}_minimal_decode_hg38.bed",
            sep="\t", header=False, index=False)
    print(f"  cleaned {popsam}: {len(df):,} segments kept")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
