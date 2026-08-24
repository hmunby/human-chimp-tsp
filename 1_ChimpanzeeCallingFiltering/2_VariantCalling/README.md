# 2_VariantCalling

GATK4 joint genotyping of the chimpanzee panel from the deduplicated BAMs (mapped to panTro6).

1. `1_haplotypecaller.sh` — `HaplotypeCaller` in GVCF mode, one gVCF per sample per chromosome
   (one array task per sample, looping over chromosomes).
2. `2_genomicsdbimport.sh` — `GenomicsDBImport` consolidates the per-sample gVCFs into a
   per-chromosome GenomicsDB workspace (one array task per chromosome). Needs a per-chromosome
   sample map (`<SET>.<chrom>.sample_map`, TSV of `sample<TAB>gvcf_path`).
3. `3_genotypegvcfs.sh` — `GenotypeGVCFs` per chromosome (heterozygosity prior 0.001) → raw
   multi-sample VCFs `genotyped.<chrom>.vcf.gz`, the input to `../3_Filtering/`.

Reference contig names are read from two lists: full reference names (for `-L`) and short names
(`1, 2A, 2B, …`) used in output filenames. Set paths/sample set in each script's parameter block.
