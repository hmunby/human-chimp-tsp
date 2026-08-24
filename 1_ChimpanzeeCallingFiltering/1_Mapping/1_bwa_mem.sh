#!/usr/bin/env bash
# Map paired-end chimpanzee reads to the panTro6 reference (GCA_002880755.3, Clint_PTRv2) with
# bwa mem, adding read groups inline, and coordinate-sort with samtools. One SLURM array task per
# sequencing run. FASTQ come from several projects (PRJEB5937, PRJEB29710, PRJDB3537,
# PRJEB15086, and samples from this study: PRJNA1357362) ; run once per project by setting PROJECT and the ID lists.
# Tools: bwa 0.7.17 (bwa mem); samtools 1.21 (sort / index / per-sample merge).

#SBATCH --job-name=bwa_mem
#SBATCH -c 8
#SBATCH --time=96:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --output logs/%x.%A_%a.out
#SBATCH --error  logs/%x.%A_%a.err
#SBATCH --array=0-8            # set to (number of runs in PROJECT) - 1
set -euo pipefail

# ---- parameters (edit) --------------------------------------------------------------------------
PROJECT="PRJEB5937"
REF="reference/GCA_002880755.3_Clint_PTRv2/GCA_002880755.3_Clint_PTRv2_genomic.fa"   # bwa-indexed
FASTQ_DIR="FASTQ/${PROJECT}/${PROJECT}_raw"
OUT_DIR="bam/mapped"
SAMPLE_IDS="FASTQ/${PROJECT}/${PROJECT}_sample_IDs.txt"   # one sample ID per line
RUN_IDS="FASTQ/${PROJECT}/${PROJECT}_run_IDs.txt"         # matching run (SRR/ERR) ID per line
# ------------------------------------------------------------------------------------------------

mkdir -p "$OUT_DIR" logs
k=$(( SLURM_ARRAY_TASK_ID + 1 ))                          # 1-based for sed
sample_ID=$(sed "${k}q;d" "$SAMPLE_IDS")
run_ID=$(sed "${k}q;d" "$RUN_IDS")
out_fn="${OUT_DIR}/${run_ID}.bam"
[ -s "$out_fn" ] && { echo "$out_fn exists; skipping"; exit 0; }

echo "Aligning run ${run_ID} (sample ${sample_ID}) -> ${out_fn}"
bwa mem "$REF" "${FASTQ_DIR}/${run_ID}_1.fastq" "${FASTQ_DIR}/${run_ID}_2.fastq" -t 8 -M \
    -R "@RG\tID:${run_ID}\tSM:${sample_ID}\tPL:ILLUMINA\tLB:${sample_ID}" \
  | samtools sort -@8 -o "$out_fn" -
echo "Done."
