# 4 — SINGER genealogical reconstruction of TSP regions

Reconstructs genealogies (ARGs) at the candidate TSP loci from stage 3 using **SINGER**, to confirm
trans-species polymorphisms and produce the final TSP set.

## Flow
1. **Group the candidate SNPs into recall windows.** The old shared SNPs (stage 3) are grouped into
   genomic windows to decide where chimpanzee variants need to be (re)called.

   > **Note — do not conflate two different "regions".** This grouping exists *only* to decide which
   > genomic windows to recall chimp variants in. It is **not** the region grouping of the final
   > TSPs we report. The final TSP regions are defined later, from the SINGER reconstruction output;
   > the recall windows here are just a coarse scaffold for variant calling.

2. **Chimp recalling** (`tree_reconstruction`): call chimpanzee variants (on hg19) within each recall
   window and merge with the human 1000GP genotypes → per-region merged VCFs (`merged_regions/`).
3. **SINGER**: downsample the human samples (keeping candidate-SNP carriers), run SINGER per
   region → ARGs → per-region TSP metrics. Run as **three chains differing only in the MCMC seed**
   (`3_SINGER`, runs `ne50000_seed1/2/3`), so the reported set can be called from pooled support rather
   than from one chain.
4. **TSP regions**: account for which candidates were testable, pool the chains into the reported
   set, cluster into regions, verify each cluster is one linked haplotype, and assign each region
   a gene.

## Components
- `1_RecallRegions/` — group candidate SNPs into recall windows.
- `2_ChimpRecall/` — GATK recall on hg19 + SHAPEIT phasing + merge with human → per-region VCFs.
- `3_SINGER/` — SINGER per region → ARGs → per-region TSP metrics. One Snakefile runs the
  three seeded chains.
- `4_TSPRegions/` — candidate-testing accounting (802/860), seed pooling to **117 TSPs**, 10 kb
  clustering to **59 regions** (**55 genes**), within-region LD verification, genetic span, and
  gene assignment.

## Outputs

Each sub-stage writes its results next to its scripts, and later stages read them from there:
`2_ClusteringLD/` produces the TSP tables and clusters, `3_GeneAssignment/` the gene
assignments, `1_CandidateTesting/` the testable-candidate accounting. Those files are bundled
as the paper run's results, so stages 5 and 6 can be run without re-running SINGER.
