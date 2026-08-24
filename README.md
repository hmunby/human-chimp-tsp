# Trans-species polymorphism between humans and chimpanzees

Analysis code for Munby & Przeworski (2026). Stages run in order; each directory has its own README with the parameters,
the inputs to set, and the software versions.

## Sub-directories

1. **`1_ChimpanzeeCallingFiltering/`** — chimpanzee variant calling and filtering: read mapping
   (bwa mem, MarkDuplicates) → GATK joint genotyping (HaplotypeCaller → GenomicsDBImport →
   GenotypeGVCFs) → filtering Snakemake → liftOver panTro6 → hg19.
2. **`2_1000GenomesProjectProcessing/`** — human variation data processing: biallelic variable sites, mappability
   masking, and annotation (allele ages, CpG status, mutation rate).
3. **`3_SharedVariation/`** — intersect human and chimpanzee variants → the shared SNP set, then
   selection of the 860 candidate TSPs (YRI midpoint Relate age ≥ 4 Mya, MHC excluded). 
4. **`4_SINGER_GenealogicalReconstruction_TSP/`** — recall chimpanzee variants at the candidate
   loci on hg19, reconstruct genealogies with SINGER, and testing of TSP status from the posterior
   ARGs. 
5. **`5_UncoalescedRegions/`** — chimpanzee demographic inference and posterior decoding (cobraa)
   → regions that have not coalesced between human and chimpanzee, their overlap,
   conservation, background-selection dependence, and gene-set enrichment.
6. **`6_TSP_Analyses/`** — characterisation of the TSP set: comparisons against the other SNP
   sets, against prior balancing-selection studies, the master table, and gene-set enrichment.

**`RESOURCES.md`** lists the external inputs each stage expects under its own `resources/`
directory — public reference data, the software that must be on `PATH`, and which stage outputs
later stages consume. Nothing in `resources/` is tracked.

`common/` holds the shared figure configuration and plotting helpers; see its README for the
`PAPER_DATA_ROOT` / `PAPER_FIG_OUTDIR` overrides.

## Conventions

- **Coordinates.** Stages 1–4 work in hg19/b37; stage 5 and the gene assignment work in hg38.
  Tables that carry both name the hg38 columns explicitly (`CHROM_HG38`, `POSITION_HG38`).
- **Ages** are stored in generations and reported in years, at 29 years per generation.
- **The TSP set** is the **pooled** set: **117** SNPs in **59** regions across **55** genes, from
  three SINGER chains combined by mean posterior support. Files carrying it are tagged `.pooled`.
  The `.strict`, `.ne20k` and `.ne100k` tags are sensitivity sets, and 860 (candidates) appears
  upstream. 
