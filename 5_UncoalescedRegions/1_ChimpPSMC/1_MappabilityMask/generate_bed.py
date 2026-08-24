#!/usr/bin/env python3
"""Convert a SNPable mask FASTA to a BED of callable regions.

SNPable encodes each position as a digit 0-3 giving how well the k-mers covering it map back
to the reference. 3 means every covering k-mer maps uniquely and without mismatch, 2 means the
majority do; 0 and 1 are unreliable. Positions scored 2 or 3 are kept as callable.

Usage: generate_bed.py <mask.fa> <out.bed>
"""
import sys

from Bio import SeqIO

CALLABLE = {"2", "3"}


def main(fasta_path, bed_path):
    with open(bed_path, "w") as bed:
        for record in SeqIO.parse(fasta_path, "fasta"):
            start = None
            for i, base in enumerate(record.seq):
                if base in CALLABLE:
                    if start is None:
                        start = i
                elif start is not None:
                    bed.write(f"{record.id}\t{start}\t{i}\n")
                    start = None
            if start is not None:   # interval running to the end of the contig
                bed.write(f"{record.id}\t{start}\t{len(record.seq)}\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
