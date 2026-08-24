# 4 / 4_TSPRegions — from confirmed TSPs to gene-assigned TSP regions

Takes the SINGER output and turns it into the reported TSP set: account for which candidates were
testable, call the TSPs by pooling the seed replicates, cluster them into regions, verify each
cluster really is one linked haplotype, and assign each region a gene.

## Flow

```
860 candidate SNPs ──1──> 802 tested for TSP status
                             │
../3_SINGER/singer/ne50000_seed{1,2,3} ┴──> 3 chains x 100 posterior trees per candidate
                                    │
                              2 ────┴──> 117 TSPs (mean support over the 3 chains > 80)
                                          └──> 59 regions (10 kb single linkage)
                                                └──> r^2 check within each multi-TSP region
                                                      │
                                                3 ────┴──> one gene per region (55 genes)
```

1. **`1_CandidateTesting/`** — the Methods accounting. **802 of 860** candidates (93.3%) were
   testable.
2. **`2_ClusteringLD/`** — pool the seed replicates into the reported set, cluster at 10 kb into
   **59 regions**, then verify LD within each of the 15 multi-TSP regions and measure their
   genetic span.
3. **`3_GeneAssignment/`** — Tiered OpenTargets molQTL / nearest-gene assignment, one gene per
   region.

## Two different things called a "region"

- The **recall region** (`../1_RecallRegions/`) is a ±250 kb window around a candidate SNP, chosen
  only to decide where to call chimpanzee variants. There are 459 of them.
- The **TSP region** (here, step 2) is the reported extent of a trans-species haplotype: the span
  of TSPs clustered at 10 kb. There are 59.

`aggregate_seeds.py` sets the `region` column to the TSP region id, keeping the recall region in
the cluster table's `orig_regions` column. Every table downstream of step 2 means TSP region.

## Coordinates

Clustering and LD are on **hg19/b37**, matching the 1000 Genomes phase 3 panel. Gene assignment
is on **hg38**, via the `CHROM_HG38`/`POSITION_HG38` columns that `../3_SINGER/` lifted. Both
coordinate systems travel in the same table; check the column name before joining.

## Scope

Table 1, the supplementary region table and the literature cross-referencing build on these outputs live in `../../6_TSP_Analyses/master_table/` and `../../6_TSP_Analyses/literature/`. Gene-set enrichment over the assigned genes, and the matched null-locus background it needs, live in `../../6_TSP_Analyses/enrichment/geneset/`.

## Software
Python 3 (numpy, pandas); bcftools 1.20; PLINK v1.90; bedtools; the OpenTargets Platform GraphQL
API.
