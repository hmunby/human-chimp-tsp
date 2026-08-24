#!/usr/bin/env bash
# Joint genotyping per chromosome from the GenomicsDB workspaces (GATK GenotypeGVCFs).
# Produces the raw multi-sample VCFs (genotyped.<chrom>.vcf.gz) that are the input to the
# filtering Snakemake in 3_Filtering/. One SLURM array task per chromosome.
# Tool: GATK 4.2.3.0.
#
#SBATCH --job-name=genotypegvcfs
#SBATCH --cpus-per-task=2
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=6G
#SBATCH --output logs/%x.%A_%a.out
#SBATCH --error  logs/%x.%A_%a.err
#SBATCH --array=0-23           # one task per chromosome
set -euo pipefail

# ---- parameters (edit) --------------------------------------------------------------------------
REF="reference/GCA_002880755.3_Clint_PTRv2/GCA_002880755.3_Clint_PTRv2_genomic.fa"
CHROM_NAMES="reference/GCA_002880755.3_Clint_PTRv2/chromnames.list"
CHROM_SIMPLE="reference/GCA_002880755.3_Clint_PTRv2/chromnames_simple.list"
SAMPLE_SET="chimp59"
GENDB_DIR="vcf/consolidated_gvcf/${SAMPLE_SET}"
OUT_DIR="vcf/genotype_calls/${SAMPLE_SET}/unfiltered"
HET="0.001"                    # -heterozygosity prior
# ------------------------------------------------------------------------------------------------

mkdir -p "$OUT_DIR" logs
readarray -t chrom  < "$CHROM_NAMES"
readarray -t cshort < "$CHROM_SIMPLE"
i=$SLURM_ARRAY_TASK_ID

gatk --java-options "-Xms4G -Xmx4G" GenotypeGVCFs \
    -R "$REF" \
    -V "gendb://${GENDB_DIR}/${cshort[$i]}" \
    -L "${chrom[$i]}" \
    -heterozygosity "$HET" \
    -O "${OUT_DIR}/genotyped.${cshort[$i]}.vcf.gz"
