#!/usr/bin/env bash
# 1000GP processing, step 1: from the 1000 Genomes phase 3 WGS *sites* VCF (b37), keep only biallelic
# SNPs.
# Tool: bcftools 1.16.
#
#SBATCH --job-name=biallelic_sites
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --output logs/%x.out
#SBATCH --error  logs/%x.err
set -euo pipefail

# ---- parameters (edit) --------------------------------------------------------------------------
# 1000GP phase3 WGS sites VCF (b37); see README for the download URL.
IN_VCF="resources/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz"
OUT_DIR="output"
# ------------------------------------------------------------------------------------------------

mkdir -p "$OUT_DIR" logs

# Keep biallelic (-m2 -M2) SNPs (-v snps).
bcftools view -m2 -M2 -v snps -Oz -o "$OUT_DIR/1000GP.wgs.biallelic.sites.vcf.gz" "$IN_VCF"
tabix "$OUT_DIR/1000GP.wgs.biallelic.sites.vcf.gz"
