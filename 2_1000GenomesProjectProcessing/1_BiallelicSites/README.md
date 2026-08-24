# 1_BiallelicSites — base 1000GP variant set

Turns the 1000 Genomes phase 3 WGS sites VCF into the biallelic, segregating SNP set that
`../3_Annotation/` annotates.

## Steps
1. `1_biallelic_sites.sh` — keep biallelic SNPs (`bcftools view -m2 -M2 -v snps`) →
   `1000GP.wgs.biallelic.sites.vcf.gz`
2. `2_variable_sites_only.sh` — drop fixed sites (`AF = 0` or `1`) →
   `1000GP.wgs.biallelic.variable_sites_only.vcf.gz`

## Inputs
- **1000GP phase 3 WGS sites VCF** (b37), release 20130502:
  `http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz`

## Software
- bcftools 1.16 (htslib 1.16); tabix.
