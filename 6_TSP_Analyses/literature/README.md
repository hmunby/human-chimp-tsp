# Literature comparisons

The published sets that the TSP regions are compared against. **The matching itself is done by
`../master_table/build_region_tables.py`**, which fills the `Literature` column of Table 1 and
Supplementary Table 2. This directory holds the reference data that it reads.

`build_region_tables.py` applies one consistent rule across all five studies.

## How each study is matched

| study | published set | matched by |
|---|---|---|
| Leffler et al. 2013 | 125 shared-polymorphism regions (hg19) | coordinate overlap |
| Rasmussen et al. 2014 | top-20 high-TMRCA regions (hg19) | coordinate overlap |
| Bitarello et al. 2018 | outlier windows only, top 0.05%, any of 12 scans (hg19) | window overlap |
| DeGiorgio et al. 2014 | top-ranked genes | TSP position inside the gene's transcribed span |
| Rocha et al. 2026 | genes deeply coalescing in **both** human and chimpanzee | as DeGiorgio |

Studies that report gene names are matched on **gene-body containment of the TSP positions** (hg38, GENCODE v49 with read-throughs stripped). The matcher records the containing gene when it differs from the assigned gene, which is where the SH3YL1 / SNTG2 footnote comes from.

## Data

| file | study | notes |
|---|---|---|
| `leffler_regions/leffler_125_regions.hg19.tsv` | Leffler 2013 | the 125-region longlist |
| `rasmussen_2014/rasmussen_2014_top20_highTMRCA_regions.hg19.tsv` | Rasmussen 2014 | |
| `bitarello_2018/bitarello_windows.hg19.tsv` | Bitarello 2018 | parsed from the supplement by `parse_bitarello_supp.py` |
| `degiorgio_2014/degiorgio_top100_genes.txt` | DeGiorgio 2014 | derived from the four pasted top-100 tables also bundled here |
| `rocha_2026/SuppTable_species_pops_overlap_summary_6MYA_NatGenetics.csv` | Rocha 2026 | the population-resolved supplementary table; the gene list is derived from it at run time, restricted to `species_sharing` of `chimp_human` or `all_three` |

