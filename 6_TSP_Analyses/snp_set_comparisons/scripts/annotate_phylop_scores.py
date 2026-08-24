"""Annotate one chromosome's SNP positions with phyloP scores.

Left-joins a CHROM/POS list onto a per-chromosome phyloP table, so positions with no score
in the alignment are kept with PHYLOP_SCORE = NA (the downstream scripts drop them). Both
inputs are hg19/b37; the phyloP table carries "chr"-prefixed names, the SNP list does not.

Run via Snakemake (uses the `snakemake` object).

Output: CHROM<TAB>POS<TAB>PHYLOP_SCORE, with a header.
"""
import pandas as pd

chrom = int(snakemake.wildcards.chrom)

snps = pd.read_csv(snakemake.input.snp_file, sep="\t", header=None, names=["CHROM", "POS"])
snps["CHROM"] = snps["CHROM"].astype(int)

phylop = pd.read_csv(snakemake.input.phylop_file, sep="\t", header=None,
                     names=["CHROM", "POS", "PHYLOP_SCORE"])
phylop["CHROM"] = phylop["CHROM"].str.replace("chr", "", regex=False).astype(int)

merged = pd.merge(snps[snps["CHROM"] == chrom], phylop, on=["CHROM", "POS"], how="left")
merged.to_csv(snakemake.output[0], sep="\t", index=False)
