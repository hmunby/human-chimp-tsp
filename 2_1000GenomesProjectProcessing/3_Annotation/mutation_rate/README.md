# Mutation rate (Roulette)

Builds `annot_table_genome.txt.gz` (`CHROM,POS,REF,ALT,MR`), used by the annotation driver to add the
`MR` INFO field. `MR` is the **Roulette** mutation-rate estimate (Seplyarskiy et al. 2023).

## Steps
1. `download_roulette.sh` — fetch the hg19 per-chromosome rate VCFs (v5.2, TFBS correction) from
   `http://genetics.bwh.harvard.edu/downloads/Vova/Roulette/hg19/autosomes/`.
2. `snakemake` (`Snakefile`) — per chrom `bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%MR\n'` →
   `annot_table_chr{N}.txt.gz`, then concatenate chr1–22 → `annot_table_genome.txt.gz` (bgzip + tabix).

hg19, autosomes only. The `MR` value is the raw Roulette rate (the source VCFs also carry `AR`
adjusted-Roulette, `MG` gnomAD, `MC` Carlson 2018, `PN` pentanucleotide — not used here).

## Caveat
`../misc_files/MutRate.hdr` declares `MR` as `Type=String` even though the values are numeric floats.
Re-declare as `Type=Float` if you need numeric filtering on `MR`.

## Software
- bcftools 1.20 (htslib 1.20); Snakemake 7.32.3; tabix.
