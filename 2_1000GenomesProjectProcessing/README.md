# 2 — 1000 Genomes Project processing

Produces the human (species-A) variant set used by `3_SharedVariation`: 1000 Genomes phase 3 SNPs,
annotated (allele ages, CpG, snpEff, mutation rate), restricted to human↔chimp mappable sites, at
MAF≥0.05 (b37), plus an hg38 liftover.

## Components
- `1_BiallelicSites/` — 1000GP phase 3 WGS sites → biallelic segregating SNPs (the base file).
- `2_Mappability/` — the human↔chimp mappable-sites list (a filter used by the annotation driver).
- `3_Annotation/` — annotate the base file → the final MAF≥0.05 species-A VCF (+ hg38 liftover).

## Data download (1000GP phase 3, release 20130502)
From the EBI FTP `http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/`:
- WGS sites VCF (v5c, b37): `ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz`
  — input to `1_BiallelicSites/`.
- Per-chromosome genotype VCFs (v5b, b37):
  `ALL.chr{1..22}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz`
  — input to `2_Mappability/`.
- Sample / population panel: `integrated_call_samples_v3.20130502.ALL.panel` (2,504 samples).

## End product
`3_Annotation/output/1000GP.wgs.biallelic.sites.YRI.CpG.snpeff.MR.MAF05.vcf.gz` (b37) — the
species-A input to `3_SharedVariation/1_Intersect/`.
