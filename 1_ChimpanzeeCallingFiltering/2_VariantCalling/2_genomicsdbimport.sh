#!/usr/bin/env bash
# Consolidate the per-sample gVCFs into a per-chromosome GenomicsDB workspace (GATK
# GenomicsDBImport). One SLURM array task per chromosome.
# Tool: GATK 4.2.3.0.
#
#SBATCH --job-name=genomicsdbimport
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=6G
#SBATCH --output logs/%x.%A_%a.out
#SBATCH --error  logs/%x.%A_%a.err
#SBATCH --array=0-23           # one task per chromosome
set -euo pipefail

# ---- parameters (edit) --------------------------------------------------------------------------
CHROM_NAMES="reference/GCA_002880755.3_Clint_PTRv2/chromnames.list"
CHROM_SIMPLE="reference/GCA_002880755.3_Clint_PTRv2/chromnames_simple.list"
SAMPLE_MAP_DIR="vcf/sample_maps"   # per-chrom sample maps: <SET>.<cshort>.sample_map (TSV: sample<TAB>gvcf_path)
SAMPLE_SET="chimp59"
OUT_DIR="vcf/consolidated_gvcf/${SAMPLE_SET}"
TMP_DIR="vcf/tmp"
# ------------------------------------------------------------------------------------------------

mkdir -p "$OUT_DIR" "$TMP_DIR" logs
readarray -t chrom  < "$CHROM_NAMES"
readarray -t cshort < "$CHROM_SIMPLE"
i=$SLURM_ARRAY_TASK_ID

gatk --java-options "-Xms4G -Xmx4G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenomicsDBImport \
    --genomicsdb-workspace-path "${OUT_DIR}/${cshort[$i]}" \
    --sample-name-map "${SAMPLE_MAP_DIR}/${SAMPLE_SET}.${cshort[$i]}.sample_map" \
    -L "${chrom[$i]}" \
    --genomicsdb-shared-posixfs-optimizations true \
    --tmp-dir "$TMP_DIR"
