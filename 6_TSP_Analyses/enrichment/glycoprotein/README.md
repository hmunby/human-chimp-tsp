# Glycoprotein gene sets

The two UniProt gene sets used by the enrichment tests, and the test for whether glycoproteins
look enriched because of where they sit in the genome rather than because of biology.

The enrichment results themselves are produced elsewhere: `../geneset/` for the TSP-assigned
genes, `../../../5_UncoalescedRegions/4_Enrichment/` for the uncoalesced regions.

## Gene sets

| set | definition | genes with coordinates |
|---|---|---|
| glycoprotein | UniProt KW-0325, human and reviewed | 4,535 |
| membrane glycoprotein | KW-0325 ∩ KW-0472 (Membrane) | 3,374 |

No single UniProt keyword means "membrane glycoprotein", so it is the intersection of the
glycoprotein and membrane keywords.

```bash
python fetch_glycoprotein_coords.py           # -> human_glycoproteins_KW0325_with_coords.tsv
python fetch_membrane_glycoprotein_coords.py  # -> human_membrane_glycoproteins_KW0325_KW0472_with_coords.tsv
```

`fetch_glycoprotein_coords.py` queries the UniProt REST API for KW-0325 and maps accessions to
hg38/hg19 coordinates through Ensembl BioMart. The membrane script subsets that output by the
KW-0325 ∩ KW-0472 accessions, so coordinates are identical to the parent set and BioMart is not
queried twice. Both `.tsv` files are bundled, so neither needs to be re-run.

## Ascertainment test

```bash
python glycoprotein_permutation_test.py       # -> glycoprotein_permutation_test.pdf, …results.tsv
```

Genes are assigned to regions by `bedtools closest`, so a gene set concentrated in gene-sparse
neighbourhoods will be picked up more often whatever the biology. This shuffles the merged
uncoalesced regions across the autosomes and recomputes the glycoprotein fraction of the nearest
genes, 1,000 times, twice: an unrestricted shuffle, and one that excludes segmental duplications.

Observed: 390 of 1,113 unique nearest genes are glycoproteins (35.0%), against a background rate
of 22.4%. Empirical p-values are on each panel of the figure.

## Not bundled

`resources/expected_uncoal/intersect/` (stage 5), the GENCODE v49 protein-coding gene bed,
`resources/hg38.genome` and `resources/chm13v2.0_SD.pos_only.hg38_coords_only.bed`.
See [`../../../RESOURCES.md`](../../../RESOURCES.md).

## Software

Python 3 (pandas, numpy, matplotlib, scipy); bedtools; network access for UniProt and BioMart.
