#!/usr/bin/env python3
"""Protein-coding gene spans with read-through transcripts excluded.

The gene-based literature criteria (DeGiorgio, Rocha) ask whether a region's TSPs fall in
the transcribed span of a listed gene. Taking that span from GENCODE's `gene` record is
wrong, because a gene record is min(start)/max(end) over ALL its transcripts -- including
read-through transcripts that continue into neighbouring genes.

    HBE1  gene record   chr11:5,268,345-5,505,652   237 kb
          canonical     chr11:5,268,345-5,269,945   1.6 kb

HBE1's two long transcripts run from the beta-globin cluster through TEN olfactory receptor
genes. Using the gene record put a TSP at chr11:5,352,508 "inside HBE1" when it is 82 kb
past the gene's real 3' end, sitting between OR51B6 and OR51M1. That single artefact put a
spurious Rocha citation on the OR51E2 region in both tables.

A gene's span is therefore taken over its transcripts EXCLUDING any that overlap another
protein-coding gene's canonical transcript. Genes with no such transcript fall back to their
canonical one.

Deliberately not a length cutoff: the genes this touches are mostly huge for real, and a
size filter cannot tell them apart. Of the eight gene-based matches in the pooled set, seven
exceed 250 kb and all seven are genuine -- RBFOX1 (2.47 Mb) really does have alternative
promoters spanning the locus, and the TSP falls in three non-read-through transcripts.

Usage:
    python build_gene_spans.py            # -> resources/gencode.v49.genes.no_readthrough.bed
"""
import argparse
import gzip
import os
import re

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
GTF = "resources/gencode.v49.annotation.gtf.gz"
OUT = os.path.join(HERE, "resources/gencode.v49.genes.protein_coding.autosomes.no_readthrough.bed")

AUTOSOMES = {f"chr{c}" for c in range(1, 23)}
_NAME = re.compile(r'gene_name "([^"]+)"')
_TYPE = re.compile(r'gene_type "([^"]+)"')


def read_transcripts(gtf):
    """Protein-coding transcripts on the autosomes: chrom, start, end, gene, canonical."""
    rows = []
    with gzip.open(gtf, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t")
            if f[2] != "transcript" or f[0] not in AUTOSOMES:
                continue
            attr = f[8]
            if (m := _TYPE.search(attr)) is None or m.group(1) != "protein_coding":
                continue
            rows.append((f[0], int(f[3]), int(f[4]), _NAME.search(attr).group(1).upper(),
                         "Ensembl_canonical" in attr))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "gene", "canon"])


def clean_spans(tx):
    """Per gene, the span over transcripts that do not run through another gene."""
    out = []
    for chrom, sub in tx.groupby("chrom", sort=False):
        canon = sub[sub["canon"]]
        cs, ce, cg = (canon["start"].values, canon["end"].values, canon["gene"].values)
        for gene, g in sub.groupby("gene", sort=False):
            keep = [r for r in g.itertuples()
                    if not ((cs < r.end) & (ce > r.start) & (cg != gene)).any()]
            if keep:
                start, end = min(r.start for r in keep), max(r.end for r in keep)
            else:                                   # every transcript trespasses
                c = g[g["canon"]]
                c = c if len(c) else g
                start, end = int(c["start"].min()), int(c["end"].max())
            out.append((chrom, start, end, gene))
    return pd.DataFrame(out, columns=["chrom", "start", "end", "gene"])


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--gtf", default=GTF)
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()

    tx = read_transcripts(a.gtf)
    print(f"  {len(tx):,} protein-coding autosomal transcripts")
    sp = clean_spans(tx).sort_values(["chrom", "start"])
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    sp.to_csv(a.out, sep="\t", header=False, index=False)
    print(f"  wrote {a.out}: {len(sp):,} genes")
    for g in ("HBE1", "RBFOX1", "SNTG2", "GPR158"):
        r = sp[sp["gene"] == g]
        if len(r):
            r = r.iloc[0]
            print(f"    {g:<8} {r.chrom}:{r.start:,}-{r.end:,}  ({r.end - r.start:,} bp)")
