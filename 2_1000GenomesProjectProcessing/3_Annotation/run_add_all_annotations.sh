#!/usr/bin/env bash
# 1000GP annotation driver. Adds background selection (BMAP), allele ages (YRI), CpG, snpEff, and
# mutation rate (MR) to the base biallelic 1000GP VCF, restricts to human<->chimp mappable sites,
# makes the MAF>=0.05 subset (the species-A input to 3_SharedVariation), and also lifts everything to hg38.
#
# Build the annotation inputs first: allele_ages/, cpg/, mutation_rate/, and ../2_Mappability/.
# Tools: bcftools 1.20 (htslib 1.20, incl. +liftover from the bcftools/score plugins); snpEff 5.1.
#
#SBATCH --job-name=add_all_annotations
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=40G
#SBATCH --output logs/%x_%A.out
#SBATCH --error  logs/%x_%A.err
set -euo pipefail

# ---- parameters (edit) --------------------------------------------------------------------------
IN_VCF="../1_BiallelicSites/output/1000GP.wgs.biallelic.variable_sites_only.vcf.gz"

# background-selection B-value map (external: CADD_bestfit B-map, Murphy et al. 2022;
# bgzipped/tabixed BED, columns CHROM,FROM,TO,BMAP)
BMAP_BED="resources/whole_genome.bmap.bed.gz"

# annotation tables (built by the sibling components)
AGE_YRI="allele_ages/allele_ages_YRI_annotation.txt.gz"
CPG_TABLE="cpg/cpg_annot_table.txt.gz"
MR_TABLE="mutation_rate/annot_table_genome.txt.gz"
MAPPABLE_SITES="../2_Mappability/output/autosomes_biallelic.sites.txt"

# INFO headers
HDR_BMAP="misc_files/BMAP.hdr"
HDR_YRI="misc_files/allele_ages_annotations_YRI.hdr"
HDR_CPG="misc_files/CpG.hdr"
HDR_MR="misc_files/MutRate.hdr"
RENAME_HG19="misc_files/rename_b37_to_hg19.txt"          # 1 -> chr1 ...

# snpEff functional annotation
SNPEFF_JAR="resources/snpEff/snpEff.jar"                 # snpEff 5.1
SNPEFF_DB="hg19"

# hg38 liftover (external references + chain)
SOURCE_REF="resources/hg19/ucsc.hg19.fasta"
TARGET_REF="resources/hg38/hg38.fa"
CHAIN="resources/chains/hg19ToHg38.over.chain"

OUT_DIR="output"
INT_DIR="intermediates"
# ------------------------------------------------------------------------------------------------
mkdir -p "$OUT_DIR" "$OUT_DIR/hg19_chroms" "$INT_DIR" tmp logs

### 1. Background selection (BMAP) ###
vcf_bmap="$INT_DIR/1000GP.wgs.biallelic.sites.BMAP.vcf.gz"
bcftools annotate -a "$BMAP_BED" -h "$HDR_BMAP" -c CHROM,FROM,TO,BMAP -Oz -o "$vcf_bmap" "$IN_VCF"
tabix "$vcf_bmap"

### 2. Allele ages (YRI) ###
vcf_YRI="$INT_DIR/1000GP.wgs.biallelic.sites.YRI.vcf.gz"
bcftools annotate -a "$AGE_YRI" -h "$HDR_YRI" -c CHROM,FROM,TO,LOWER_AGE_YRI,UPPER_AGE_YRI,AGE_ESTIMATE_YRI -Oz -o "$vcf_YRI" "$vcf_bmap"
tabix "$vcf_YRI"

### 3. CpG ###
vcf_CpG="$INT_DIR/1000GP.wgs.biallelic.sites.YRI.CpG.vcf.gz"
bcftools annotate -a "$CPG_TABLE" -h "$HDR_CPG" -c CHROM,POS,REF,ALT,CPG -Oz -o "$vcf_CpG" "$vcf_YRI"
tabix "$vcf_CpG"

### 4. snpEff (functional / loss-of-function) ###
vcf_snpeff="$INT_DIR/1000GP.wgs.biallelic.sites.YRI.CpG.snpeff.vcf"
java -Xmx4g -jar "$SNPEFF_JAR" "$SNPEFF_DB" "$vcf_CpG" > "$vcf_snpeff"
bgzip "$vcf_snpeff"
tabix "${vcf_snpeff}.gz"

### 5. Mutation rate (MR) ###
vcf_mr="$INT_DIR/1000GP.wgs.biallelic.sites.YRI.CpG.snpeff.MR.vcf.gz"
bcftools annotate -a "$MR_TABLE" -h "$HDR_MR" -c CHROM,POS,REF,ALT,MR -Oz -o "$vcf_mr" "${vcf_snpeff}.gz"
tabix "$vcf_mr"

### 6. Mappability filter (human<->chimp mappable sites) ###
vcf_mappable="$OUT_DIR/1000GP.wgs.biallelic.sites.YRI.CpG.snpeff.MR.vcf.gz"
bcftools view -R "$MAPPABLE_SITES" -Oz -o "$vcf_mappable" "$vcf_mr"
tabix "$vcf_mappable"

### 7. MAF subsets (b37). The MAF>=0.05 file is the species-A input to 3_SharedVariation. ###
bcftools view -i 'AF>=0.05 && AF<=0.95' -Oz -o "${vcf_mappable%.vcf.gz}.MAF05.vcf.gz" "$vcf_mappable"
tabix "${vcf_mappable%.vcf.gz}.MAF05.vcf.gz"
bcftools view -i 'AF>=0.04 && AF<=0.96' -Oz -o "${vcf_mappable%.vcf.gz}.MAF04.vcf.gz" "$vcf_mappable"
tabix "${vcf_mappable%.vcf.gz}.MAF04.vcf.gz"

### 8. Liftover to hg38 ###
# Rename b37 (1) -> hg19 (chr1); store hg19 coords in the ID so they survive liftover.
vcf_hg19="$OUT_DIR/hg19_chroms/1000GP.wgs.biallelic.sites.YRI.CpG.snpeff.MR.vcf.gz"
bcftools annotate --rename-chrs "$RENAME_HG19" "$vcf_mappable" -Oz -o "$vcf_hg19"
tabix "$vcf_hg19"
bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' "$vcf_hg19" -Oz -o "${vcf_hg19%.vcf.gz}.withID.vcf.gz"
tabix "${vcf_hg19%.vcf.gz}.withID.vcf.gz"

vcf_hg38="$OUT_DIR/1000GP.wgs.biallelic.sites.YRI.CpG.snpeff.MR.hg38.vcf.gz"
reject="$OUT_DIR/hg19_chroms/1000GP.wgs.biallelic.sites.YRI.CpG.snpeff.MR.LIFT_REJECTED.vcf.gz"
bcftools +liftover -Ou "${vcf_hg19%.vcf.gz}.withID.vcf.gz" -- \
    -s "$SOURCE_REF" -f "$TARGET_REF" -c "$CHAIN" \
    -Oz --reject "$reject" | bcftools sort --temp-dir tmp/ -Oz -o "$vcf_hg38"
tabix "$vcf_hg38"
tabix "$reject"

# hg38 MAF subsets
bcftools view -i 'AF>=0.05 && AF<=0.95' "$vcf_hg38" -Oz -o "${vcf_hg38%.vcf.gz}.MAF05.vcf.gz"
tabix "${vcf_hg38%.vcf.gz}.MAF05.vcf.gz"
bcftools view -i 'AF>=0.04 && AF<=0.96' "$vcf_hg38" -Oz -o "${vcf_hg38%.vcf.gz}.MAF04.vcf.gz"
tabix "${vcf_hg38%.vcf.gz}.MAF04.vcf.gz"
