# 3_SharedVariation / 1_Intersect

Intersects the human and chimpanzee SNP sets to produce the shared-variation call set, then removes
**human** paralogous regionz.

- **Species A** = human 1000GP annotated biallelic sites (b37) — output of `2_1000GenomesProjectProcessing/`.
- **Species B** = chimpanzee calls lifted to the human reference (b37) — output of
  `1_ChimpanzeeCallingFiltering/3_Filtering/` (`8_human_ref/genotyped.chimp43.hg19.b37.vcf.gz`).

## Flow (rules, in order)
1. `intersect` — `bcftools isec -n=2 -w 1,2`: sites present in both VCFs, keeping each species' records.
2. `update_extract_AF` — recompute chimp AF and extract a CHROM/POS/AF/ID annotation table.
3. `annoate_AF_ID` — annotate the human shared records with chimp AF + panTro6 ID → `shared_variants/`.
4. `filter_paralogs_SD`, `filter_paralogs_CNV_DUP` — drop sites in human SD and human CNV/DUP regions.
5. `output_variant_table` — final table (AF, chimp AF/ID, ancestral allele, CpG, YRI Relate ages, ANN).

## Files
- `Snakefile`, `config.yaml` — the pipeline and its SLURM (cluster) profile; profile is cluster-specific.
- `misc_files/` — the chimp AF/ID header + the two human paralog BEDs (see its README).

## Inputs to set (top of `Snakefile`)
- `SPECIES_A_VCF` — human annotated sites from stage 2.
- `SPECIES_B_VCF` — chimp lifted VCF from stage 1 (relative path pre-filled).

## Software
- Snakemake 7
- bcftools 1.20 (htslib 1.20) — `isec`, `+fill-tags`, `annotate`, `query`, `stats`, `view`
- tabix / bgzip (htslib 1.20)
