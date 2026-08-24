# 6 — TSP analyses

Downstream characterisation of the TSP set from stage 4: how it compares to the SNP sets it is
drawn from, how it compares to prior balancing-selection studies, the per-SNP master table and the
manuscript region tables, and gene-set enrichment over the assigned genes.

## Components

- `snp_set_comparisons/` — the TSPs against the SNP sets they are drawn from (1000GP MAF≥5%,
  all shared, old shared ≥4 Mya), plus chimpanzee for conservation: overlap with the uncoalesced
  regions (**Figure 2**), snpEff variant-effect distribution (**Panel A**), phyloP eCDF
  (**Panel B**) and the accompanying Mann-Whitney / KS tests. See its README.
- `literature/` — overlap of the TSP regions with prior studies (Leffler 2013, DeGiorgio 2014,
  Rasmussen 2014, Bitarello 2018, Soni 2022, Rocha 2026). See its README. Gene-based studies are
  matched by gene-body overlap; region-based studies by coordinate overlap. The matching that
  feeds the manuscript tables lives in `master_table/build_region_tables.py`
- `master_table/` — the manuscript tables.
  - `build_tsp_master_table.py` → `tsp_master_table.tsv`, one row per TSP, joining the TSP set
    with gene assignments, molQTL evidence and phyloP.
  - `build_region_tables.py` → `Table1_multiTSP.pooled.*` and `SuppTable_all_regions.pooled.*`,
    one row per region, including the literature column.
  - `build_gene_spans.py` → the gene-span reference used by the region tables.
- `examples/` — worked single-locus examples.
  - `muc7_locus_figure/` — **Figure 5**: the local genealogy at two of the MUC7 TSPs
    (`plot_panel_a.py`) and the per-SNP TSP-support track for the LocusZoom panel
    (`build_panel_b_input.py`). See its README.
  - `muc7_eqtl/` — eQTL evidence for the same locus. See its README.
  - `muc7_archaic/` — LD between the MUC7 TSPs and the Xu et al. (2017) archaic haplotype.
- `enrichment/` — gene-set enrichment analyses.
  - `geneset/` — GO / KEGG / Reactome / Hallmark enrichment of the assigned genes against an
    ascertainment-matched null (`build_background.py` builds the matched null loci,
    `enrichment_scan.py` runs the scan), plus the glycoprotein tests.
  - `glycoprotein/` — the glycoprotein (UniProt KW-0325) and membrane-glycoprotein
    (KW-0325 ∩ KW-0472) gene sets, and the ascertainment permutation test. See its README.

## Inputs

These scripts read the TSP set and the region/gene assignments (stage 4 outputs) and the shared-SNP
set (stage 3), plus the bundled per-study literature data. Files carrying the reported set are
tagged `.pooled`; `.strict`, `.ne20k` and `.ne100k` are the sensitivity sets. Each script has an
editable block at the top for the paths it needs.

