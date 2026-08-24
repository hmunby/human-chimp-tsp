#!/usr/bin/env bash
# 1000GP processing, step 2: drop monomorphic sites (AF = 0 or 1) to keep only segregating biallelic
# SNPs. A small number of 1000GP sites are fixed (AF 0/1), e.g. after related-sample removal upstream.
# The result is the base file that 3_Annotation/ annotates.
# Tool: bcftools 1.16.
#
#SBATCH --job-name=variable_sites_only
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH --output logs/%x.out
#SBATCH --error  logs/%x.err
set -euo pipefail

OUT_DIR="output"
IN="$OUT_DIR/1000GP.wgs.biallelic.sites.vcf.gz"
OUT="$OUT_DIR/1000GP.wgs.biallelic.variable_sites_only.vcf.gz"

bcftools view -i 'INFO/AF != 0 && INFO/AF != 1' -Oz -o "$OUT" "$IN"
tabix "$OUT"
bcftools stats "$OUT" > "${OUT%.vcf.gz}.stats"
