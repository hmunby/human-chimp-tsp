#!/usr/bin/env bash
# GATK HaplotypeCaller in GVCF mode: one gVCF per sample per chromosome, from the deduplicated
# BAMs. One SLURM array task per sample; loops over chromosomes internally.
# Tool: GATK 4.2.3.0.
#
#SBATCH --job-name=haplotypecaller
#SBATCH --cpus-per-task=2
#SBATCH --time=96:00:00
#SBATCH --mem-per-cpu=24G
#SBATCH --output logs/%x.%A_%a.out
#SBATCH --error  logs/%x.%A_%a.err
#SBATCH --array=0-21           # set to (number of samples) - 1
set -euo pipefail

# ---- parameters (edit) --------------------------------------------------------------------------
REF="reference/GCA_002880755.3_Clint_PTRv2/GCA_002880755.3_Clint_PTRv2_genomic.fa"
IN_DIR="bam/deduped"
OUT_DIR="vcf/gvcf"
SAMPLE_IDS="bam/sample_IDs.txt"
CHROM_NAMES="reference/GCA_002880755.3_Clint_PTRv2/chromnames.list"          # reference contig names (-L)
CHROM_SIMPLE="reference/GCA_002880755.3_Clint_PTRv2/chromnames_simple.list"  # short names for filenames (1,2A,...)
# ------------------------------------------------------------------------------------------------

mkdir -p "$OUT_DIR" logs
readarray -t sample < "$SAMPLE_IDS"
readarray -t chrom  < "$CHROM_NAMES"
readarray -t cshort < "$CHROM_SIMPLE"
s=${sample[$SLURM_ARRAY_TASK_ID]}

for i in "${!chrom[@]}"; do
    gatk --java-options "-Xms20G -Xmx20G -XX:ParallelGCThreads=2" HaplotypeCaller \
        -R "$REF" \
        -I "${IN_DIR}/${s}.bam" \
        -L "${chrom[$i]}" \
        -ERC GVCF \
        -O "${OUT_DIR}/${s}.${cshort[$i]}.g.vcf.gz"
done
