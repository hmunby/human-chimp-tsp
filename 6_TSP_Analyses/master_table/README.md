# Master table of the trans-species-polymorphism (TSP) SNPs

`tsp_master_table.tsv` — one row per TSP SNP (117 rows for the reported `pooled` set), 31 columns,
tab-separated. Built by `build_tsp_master_table.py`, which **joins pre-existing per-SNP tables**;
nothing is recomputed. Reproduce with:

```
python build_tsp_master_table.py     # -> tsp_master_table.tsv
```

Pass `--tsps`/`--gene-assignments`/`--chains 1` to build the same table for a sensitivity set
(`strict`, `ne20k`, `ne100k`).

---

## The TSP set

The 117 SNPs are the manuscript TSP set: shared human–chimp SNPs supported as a first- or
second-coalescence trans-species configuration in more than 80 % of the posterior SINGER
genealogies, pooled over three seed-replicate chains, then clustered to 10 kb into 59 regions.
They come verbatim from `../../4_SINGER_GenealogicalReconstruction_TSP/4_TSPRegions/2_ClusteringLD/tsps.pooled.hg19.txt`,
which already carries most fields (hg19+hg38 coordinates, chimp panTro6 id, allele frequencies,
CpG, BMAP, Roulette mutation rate, YRI age bounds in generations, and the ARG-support counts).

## Sources and join keys

| source | join key | fields taken |
|---|---|---|
| `…/2_ClusteringLD/tsps.pooled.hg19.txt` | (base table, 117 rows) | coords hg19/hg38, `CHIMP_ID`, `REF`/`ALT`/`AA`, `AF`, `CHIMP_AF`, `CPG`, `BMAP`, `MR`, `LOWER/UPPER_AGE_YRI`, `Annotation`, `tsp_count`, `tsp_first_or_second_count`, `region` |
| `…/3_GeneAssignment/variant_gene_assignments.pooled.tsv` | hg38 `CHROM_HG38`+`POSITION_HG38` | `rsid`, `most_severe`, `nearest_gene`, `nearest_distance`, `assigned_gene`, `assignment_tier`, `in_ot_index`, `n_credible_sets`, `molqtl_details`, `gwas_l2g_genes` |
| `phylop/` (36-mammal EPO, human excluded) | hg19 `CHROM`+`POS` | **phyloP, human-excluded (main-figure source)** |

All 117 TSPs match the gene-assignment table on hg38 coordinate (1-to-1).

## phyloP source (matches the main figure)

One phyloP column, the **same track used in the main phyloP figure and the K-S tests**:

- **`phyloP (36-mammal EPO, human excluded)`** — 36-eutherian-mammal EPO alignment with the human
  sequence excluded (so the score reflects conservation in the *other* mammals, chimp retained).
  This is the `phylop/` directory, the main-figure / `phylop_stats.py` source. **2 of 117** SNPs
  have no score here and are written as **`NA`**, because those positions fall outside
  EPO-alignment coverage.

## Units and conventions

- **Ages** are stored in the source in **generations** and reported here in **years**
  (years = generations × 29, the manuscript generation time). **Midbranch = (lower + upper) / 2**.
  All 117 have ages.
- **ARG support for TSP (% of 300 posterior genealogies, 3 chains × 100)** = `mean_support`, the
  mean of the three chains' `tsp_first_or_second_count`. Because each chain contributes exactly
  100 genealogies, that mean is the percentage of all 300.
- **Chimp coordinates** are panTro6, parsed from `CHIMP_ID` (`chrom_pos_ref_alt`); the ref/alt in
  that id are the chimp alleles.
- **Ancestral allele** is the 1000 Genomes ancestral call, upper-cased. (A lowercase call in the
  raw source denotes a low-confidence assignment; that flag is not carried into this table.)
- **molQTL type(s) / tissue(s)** are the unique, order-preserving set of QTL types and tissues
  across *all* colocalising credible sets for the variant (Open Targets `molqtl_details`), not just
  the single best one. Multiple values are `;`-separated (tissue names can themselves contain commas).
- **CpG status**: `CpG` / `non-CpG` from the source `CPG` 1/0 flag.
- **BMAP** is the McVicker B value (higher = weaker background selection).
- Coordinates are 1-based. The set is already MHC-excluded upstream.

## Sparse-by-nature columns (not missing data)

- `rsID` and `Most severe consequence (VEP)`: 75/117 — the other 42 are novel/un-rsID'd variants.
- `molQTL type(s)` / `molQTL tissue(s)`: 53/117 have a colocalising molQTL.
- `GWAS L2G effector genes`: 32/117 have a GWAS L2G link.

## Columns (in order)

`rsID`, `Chromosome (hg19)`, `Position (hg19)`, `Chromosome (hg38)`, `Position (hg38)`,
`Chimp chromosome (panTro6)`, `Chimp position (panTro6)`, `Reference allele`, `Alternate allele`,
`Ancestral allele`, `Ancestral low-confidence`, `Allele frequency (1000 Genomes)`,
`Allele frequency (chimpanzee)`, `CpG status`, `Background-selection B value (BMAP)`,
`Mutation rate (Roulette)`, `phyloP (36-mammal EPO, human excluded)`,
`phyloP (46-way, human included)`, `Variant annotation (snpEff)`, `Most severe consequence (VEP)`,
`Nearest gene`, `Distance to nearest gene (bp)`, `Assigned gene`, `Assignment tier`,
`In Open Targets index`, `Number of credible sets`, `molQTL type(s)`, `molQTL tissue(s)`,
`GWAS L2G effector genes`, `Lower age bound, YRI (years)`, `Upper age bound, YRI (years)`,
`Midbranch age, YRI (years)`, `ARG support for TSP (% of 300 posterior genealogies, 3 chains x 100)`.

## Literature column

Matching criteria differ by study, following the manuscript:

| study | criterion | source |
|---|---|---|
| Leffler 2013, Rasmussen 2014 | coordinate overlap with a reported interval (hg19) | `literature/*/` |
| Bitarello 2018 | overlap with an **outlier** window only, not the wider candidate set (hg19) | `literature/bitarello_2018/` |
| DeGiorgio 2014 | region's TSPs fall inside a listed gene's transcribed span (hg38, GENCODE v49) | top-100 gene lists |
| Rocha 2026 | as DeGiorgio, against the published 6-Mya supplementary table | see below |

**Rocha.** The gene list is derived by `rocha_genes()` from
`literature/rocha_2026/SuppTable_species_pops_overlap_summary_6MYA_NatGenetics.csv`
(sudmantlab/panpan_diversity_project), restricted to genes reported as deeply coalescing in
**both human and chimpanzee** — `species_sharing` of `chimp_human` or `all_three`, 29 genes.

