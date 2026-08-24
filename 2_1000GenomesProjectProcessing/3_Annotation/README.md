# 3_Annotation — annotate the 1000GP base VCF

Adds allele ages, CpG, functional (snpEff) and mutation-rate annotations to the base biallelic 1000GP
VCF (`../1_BiallelicSites/`), restricts to human↔chimp mappable sites (`../2_Mappability/`), and
produces the MAF≥0.05 species-A VCF used by `3_SharedVariation` — plus an hg38 liftover.

Driver: `run_add_all_annotations.sh`. Order: BMAP → YRI ages → CpG → snpEff → MR →
mappability filter → MAF subsets → hg38 liftover.

## Annotation inputs (build these first)
| INFO field(s) | component | source |
|---|---|---|
| `BMAP` | driver (external BED) | background-selection B-map, CADD_bestfit (Murphy et al. 2022) |
| `LOWER/UPPER/AGE_ESTIMATE_YRI` | `allele_ages/` | Relate ages (YRI), Speidel 2019 (Zenodo 3234689) |
| `CPG` | `cpg/` | ancestral-corrected CpG (MutationalPatterns + hs37d5) |
| `ANN` | snpEff (in the driver) | snpEff 5.1, hg19 db |
| `MR` | `mutation_rate/` | Roulette v5.2 TFBS (Seplyarskiy 2023) |
| — (filter, not a field) | `../2_Mappability/` | hg19↔panTro6 mappable sites |

## Output
- `output/1000GP.wgs.biallelic.sites.YRI.CpG.snpeff.MR.MAF05.vcf.gz` (b37) — the species-A input
  to `3_SharedVariation/1_Intersect/`.
- `output/…hg38.vcf.gz` (+ MAF04/MAF05 subsets) — hg38 liftover for downstream (e.g. SINGER).

## Files
- `run_add_all_annotations.sh` — the driver.
- `allele_ages/`, `cpg/`, `mutation_rate/` — build the annotation tables
- `misc_files/` — INFO headers (BMAP/YRI/CpG/MR) + the b37→hg19 chromosome rename table.

## Inputs 
- hg19 + hg38 FASTAs; the `hg19ToHg38.over.chain` liftOver chain.

## Software
- bcftools 1.20 (htslib 1.20, incl. `+liftover`); snpEff 5.1; R 4.5.1 (CpG); Snakemake 7; tabix/bgzip.
