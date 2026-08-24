# 4 / 2_ChimpRecall — recall chimp variants in the candidate regions & merge with human

The chimp recalling pipeline. For each recall window (from `../1_RecallRegions/`), it calls chimp
variants on hg19, phases them, and merges with the human 1000GP genotypes → per-region VCFs that are
the input to SINGER (`../3_SINGER/`, `merged_regions/{chrom}/{region}.vcf.gz`).

## Flow (rules, in order)
- **remap** / **index**: remap the panTro6 chimp BAMs to hg19 (`samtools fastq | bwa mem`), then
  index. This is the chimp→hg19 mapping, distinct from the chimp→panTro6 mapping of stage 1.
- **haplotype_caller_region → genomicsdb_import → genotype_gvcfs**: GATK joint genotyping per chrom
  over the recall intervals (±250 kb interval padding), with the chimp pedigree.
- **biallelic_filter → filter → setIDs**: biallelic SNPs, site/genotype QC via `vcfCleaner gatk`
  (same settings as stage 1), assign `CHROM_POS_REF_ALT` IDs.
- **compute_thin_exclude → vcf2plink → bed2ped → run_shapeit → correct_fixed → convert_hap_vcf →
  transfer_ann**: SHAPEIT phasing (duoHMM; uniform rho, no genetic map), protecting shared-SNP sites
  from thinning.
- **remove_related_chimps**: subset to the unrelated chimps.
- **subset_human → merge_vcfs → chimp_fail → filter_fill_merged → split_regions**: subset 1000GP to
  the same windows, merge human + chimp, fill missing GTs, split into per-region VCFs → `merged_regions/`.

## Inputs to set (top of `Snakefile`, `resources/`)
- `REF` — human hg19 FASTA; `BAM_DIR` — the source panTro6 chimp BAMs that `remap` reads
  (it writes the hg19 ones to `bam/`); `TGP_VCF_DIR` — 1000GP phase3
  per-chrom genotype VCFs (v5b); `SHARED_SITES` — shared-SNP sites from stage 3 (relative path pre-filled).
- Recall-window region files are bundled: `misc_files/regions_chr*.list` (GATK intervals) and
  `region_files/regions_chr*_250kb*.txt` (merged windows). See `../1_RecallRegions/`.

## misc_files 
Chimp sample→BAM map, pedigree `.fam`, unrelated-chimp list, candidate SNP positions
(`merged_candidate_snps.chrom_pos.txt`), and the per-chrom recall region lists.

## Software
- GATK 4 (as stage 1, 4.2.3.0); vcfCleaner (tplinderoth/ngsQC, `gatk` mode); PLINK v1.90b7;
  SHAPEIT v2 (duoHMM); bcftools 1.20 (htslib 1.20); bwa 0.7.17 + samtools (remap); Snakemake 7.
