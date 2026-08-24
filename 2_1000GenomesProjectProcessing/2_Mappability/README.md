# 2_Mappability — human↔chimp mappable sites

Snakemake pipeline that restricts 1000GP autosomal variants to positions that map cleanly between
human (hg19) and chimp (panTro6). The output site list is used by `../3_Annotation/` to filter the
annotated VCF.

## Criterion
A variant is kept iff it passes all three:
1. the **site** lifts hg19 → panTro6 (`bcftools +liftover`, not rejected);
2. a **±100 bp window** around it lifts contiguously (UCSC `liftOver` on a `bedtools slop` BED);
3. it is a **biallelic SNP**.

Surviving variants are reported on the original hg19 genome. This mirrors the region-lift criterion
applied to chimp variants in stage 1 (`1_ChimpanzeeCallingFiltering/3_Filtering`, rule `region_lift`).

## Files
- `Snakefile`, `config.yaml` (SLURM profile — cluster-specific).
- `misc_files/chr_b37_hg19.txt` — `1 → chr1` chromosome rename map.

## Inputs 
- 1000GP phase3 per-chromosome genotype VCFs (v5b, b37).
- hg19 + panTro6 FASTAs, hg19 chrom sizes, the `hg19ToPanTro6.over.chain` liftOver chain.

## Software
- bcftools 1.16 (with the `+liftover` bcftools/score plugin); UCSC `liftOver`; bedtools v2.29.1;
  Snakemake 7; tabix.