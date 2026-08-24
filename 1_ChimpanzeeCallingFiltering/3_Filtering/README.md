# 3_Filtering — chimpanzee VCF processing

Snakemake pipeline that turns the raw joint-genotyped chimp VCFs (`genotyped.<chrom>.vcf.gz` from
`../2_VariantCalling/`) into the filtered, hg19-lifted chimp variant set used downstream. 

## Stages (rules, in order)
1. Biallelic SNPs filter
2. Site / genotype QC filters
3. Mendelian-consistency filter
4. Subset to unrelated founders
5. Segmental-duplication filter (chimp T2T SD track)
6. Allele-count filter
7. Combine chromosomes
8. LiftOver panTro6 → hg19

## Files
- `Snakefile` — the pipeline.
- `config.yaml` — Snakemake **cluster (SLURM) profile**; cluster-specific, edit for your scheduler.
- `misc_files/` — the small non-downloadable inputs the Snakefile needs:
  - `chimp59_pedigree.fam`, `chimp59_founder_sample_IDs.txt` — pedigree + unrelated-founder list.
  - `chr_chimp.txt`, `chr_hg19_b37.txt`, `chr_b37_hg19.txt`, `hg19_chrom_names.txt` — chromosome-name maps.
  - `panTro6.chrom_lengths` — contig lengths.
- `segmental_duplication_track/` — builds the stage-5 segmental-duplication mask (BISER: T2T vs panTro6 assemblies).

## Inputs to set
Edit the `#### INPUT PATHS ####` block at the top of the `Snakefile`:
- `RAW_VCF_DIR` — raw joint-genotyped VCFs from `../2_VariantCalling/` (`genotyped.chr<chrom>.vcf.gz`).
- `resources/` — external data: panTro6 + hg19 FASTAs and the UCSC panTro6↔hg19 liftOver chains.

## Software
- Snakemake 7
- bcftools 1.20 (htslib 1.20); `+liftover` from the bcftools/score plugins (https://github.com/freeseek/score)
- [vcfCleaner](https://github.com/tplinderoth/ngsQC/tree/master/vcfCleaner)
- PLINK v1.90b7 (Mendelian-error analysis)
- tabix / bgzip (htslib 1.20)