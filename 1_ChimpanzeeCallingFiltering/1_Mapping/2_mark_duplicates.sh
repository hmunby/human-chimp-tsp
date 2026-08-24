#!/usr/bin/env bash
# Mark PCR/optical duplicates per sample with GATK/Picard MarkDuplicates (GATK 4.2.3.0).
# Input: one coordinate-sorted BAM per sample (per-run BAMs from 1_bwa_mem.sh merged per sample
# with `samtools merge` where a sample has >1 run; read groups already set inline by bwa mem).
# Output: the deduplicated per-sample BAMs that feed variant calling (2_VariantCalling/).
#
#SBATCH --job-name=mark_dups
#SBATCH -c 1
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --output logs/%x.%A_%a.out
#SBATCH --error  logs/%x.%A_%a.err
#SBATCH --array=0-9            # set to (number of samples) - 1
set -euo pipefail

# ---- parameters (edit) --------------------------------------------------------------------------
IN_DIR="bam/merged"           # one BAM per sample (post bwa mem + per-sample merge)
OUT_DIR="bam/deduped"
TMP_DIR="bam/tmp"
SAMPLE_IDS="bam/sample_IDs.txt"     # one sample ID per line
# ------------------------------------------------------------------------------------------------

mkdir -p "$OUT_DIR" "$TMP_DIR" logs
readarray -t sample < "$SAMPLE_IDS"
s=${sample[$SLURM_ARRAY_TASK_ID]}

gatk --java-options "-Xmx8G" MarkDuplicates \
    -I "${IN_DIR}/${s}.bam" \
    -O "${OUT_DIR}/${s}.bam" \
    -M "${OUT_DIR}/${s}.marked_dup_metrics.txt" \
    -MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 60000 \
    --TMP_DIR "$TMP_DIR"
samtools index "${OUT_DIR}/${s}.bam"
