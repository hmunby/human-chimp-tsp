#!/usr/bin/env bash
# Segmental-duplication track, step 1: run BISER on the (soft-masked) T2T chimpanzee assembly
# (mPanTro3 hap1) and the panTro6 reference together.
# Tool: BISER v1.4  (pip install biser; https://github.com/0xTCG/biser).
#
#SBATCH --job-name=biser
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --output logs/%x.%A.out
#SBATCH --error  logs/%x.%A.err
set -euo pipefail

# ---- parameters (edit) --------------------------------------------------------------------------
# Soft-masked FASTAs (see README, step 0). panTro6: UCSC soft-masked download. T2T:
# mPanTro3.hap1.cur.20231122 soft-masked with the assembly's RepeatMasker annotation
# (mPanTro3_v2.0.RepeatMasker_v1.1.hap1) via `bedtools maskfasta`.
T2T="reference/T2T/mPanTro3.hap1.cur.20231122.masked.fasta"
PANTRO6="reference/GCA_002880755.3_Clint_PTRv2/UCSC/panTro6.fa"
OUT="out/t2t_pantro6_seg_dups.txt"
THREADS=4
# ------------------------------------------------------------------------------------------------

mkdir -p out logs
# --keep-contigs: panTro6 sequences are named as scaffolds/contigs, not "chrN"; without it BISER
# reports "No chromosomes found" and skips them. --gc-heap caps the SD-decomposition heap.
biser -o "$OUT" -t "$THREADS" --gc-heap 12G --keep-contigs --keep-temp "$T2T" "$PANTRO6"
