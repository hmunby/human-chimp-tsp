#!/usr/bin/env python3
"""Merge bookended (adjacent) BED entries and compute length-weighted averages of score columns.

Adjacent entries are those where end[i] == start[i+1] on the same chromosome.
All score columns (col 4 onward) are averaged weighted by interval length.

Snakemake usage:
    rule merge_bed:
        input:  bed="{prefix}.bed"
        output: bed="{prefix}.merged.bed"
        script: "merge_bed_weighted.py"

Standalone usage:
    python merge_bed_weighted.py input.bed -o output.bed
"""

import sys
import argparse


def parse_scores(raw):
    """Parse score fields; NA stays as None."""
    return [None if x == "NA" else float(x) for x in raw]


def merge_bed_weighted(input_file, output_file):
    current_chrom = None
    current_start = None
    current_end = None
    weighted_sums = None   # per-column sum of (value * length), None if all-NA so far
    col_weights = None     # per-column total length of non-NA entries
    total_length = 0

    def flush(out):
        cols = [current_chrom, str(current_start), str(current_end)]
        for ws, cw in zip(weighted_sums, col_weights):
            cols.append("NA" if cw == 0 else str(ws / cw))
        out.write("\t".join(cols) + "\n")

    def init(chrom, start, end, scores, length):
        ws = [0.0 if s is None else s * length for s in scores]
        cw = [0 if s is None else length for s in scores]
        return chrom, start, end, ws, cw, length

    with open(input_file) as f_in:
        raw_lines = [l.rstrip("\n") for l in f_in if l.strip() and not l.startswith("#")]

    def sort_key(l):
        p = l.split("\t")
        chrom = p[0]
        num = int(chrom[3:]) if chrom[3:].isdigit() else float("inf")
        return (num, chrom, int(p[1]))

    raw_lines.sort(key=sort_key)

    with open(output_file, "w") as f_out:
        for line in raw_lines:
            parts = line.split("\t")
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            scores = parse_scores(parts[3:])
            length = end - start

            if current_chrom is None:
                current_chrom, current_start, current_end, weighted_sums, col_weights, total_length = init(
                    chrom, start, end, scores, length
                )
            elif chrom == current_chrom and start == current_end:
                current_end = end
                total_length += length
                for i, s in enumerate(scores):
                    if s is not None:
                        weighted_sums[i] += s * length
                        col_weights[i] += length
            else:
                flush(f_out)
                current_chrom, current_start, current_end, weighted_sums, col_weights, total_length = init(
                    chrom, start, end, scores, length
                )

        if current_chrom is not None:
            flush(f_out)


try:
    input_file = snakemake.input.bed
    output_file = snakemake.output.bed
    merge_bed_weighted(input_file, output_file)
except NameError:
    parser = argparse.ArgumentParser(
        description="Merge adjacent BED entries and compute length-weighted averages of score columns."
    )
    parser.add_argument("input", help="Input BED file")
    parser.add_argument("-o", "--output", required=True, help="Output BED file")
    args = parser.parse_args()
    merge_bed_weighted(args.input, args.output)
