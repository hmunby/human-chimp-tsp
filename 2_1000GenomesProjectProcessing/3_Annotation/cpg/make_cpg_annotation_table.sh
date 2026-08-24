#!/usr/bin/env bash
# Build the per-chromosome CpG annotation tables (one SLURM array task per autosome), then combine
# with combine_cpg_tables.sh. Each table is CHROM,POS,REF,ALT,CPG where CPG=1 marks an ancestral
# C->T mutation in a CpG dinucleotide context (see README + make_cpg_annotation.R).
#
#SBATCH --job-name=make_CpG_annot
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --output logs/%x_%A_%a.out
#SBATCH --error  logs/%x_%A_%a.err
#SBATCH --array=1-22
set -euo pipefail

# ---- parameters (edit) --------------------------------------------------------------------------
# Base biallelic 1000GP VCF from ../../1_BiallelicSites/.
VCF_GENOME="../../1_BiallelicSites/output/1000GP.wgs.biallelic.variable_sites_only.vcf.gz"
CONDA_ENV="mut_patterns"     # env with R 4.5.1 + MutationalPatterns/VariantAnnotation/BSgenome.hs37d5
# ------------------------------------------------------------------------------------------------

i=$SLURM_ARRAY_TASK_ID
mkdir -p tmp cpg_tables_per_chrom logs

# Per-chromosome VCF
bcftools view -r "${i}" "$VCF_GENOME" -Oz -o "tmp/chr${i}.vcf.gz"
tabix "tmp/chr${i}.vcf.gz"

# CpG annotation table for this chromosome
mamba activate "$CONDA_ENV"
Rscript make_cpg_annotation.R "tmp/chr${i}.vcf.gz" "cpg_tables_per_chrom/cpg.annot_table.chr${i}.txt"

rm -f "tmp/chr${i}.vcf.gz" "tmp/chr${i}.vcf.gz.tbi"
