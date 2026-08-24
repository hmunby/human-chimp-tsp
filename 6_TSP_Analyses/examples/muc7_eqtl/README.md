# MUC7 eQTL evidence

The five TSPs of the MUC7 region lie in the intergenic sequence between *CABS1* and *SMR3A*, but
fall in the credible set of an eQTL for *MUC7* in skin. Minor salivary gland is the only salivary
tissue with eQTL data, and there the TSPs show no association with *MUC7* expression despite the
gene being highly expressed. These scripts are where that comes from.

## The variants

| GRCh38 | GRCh37 | rsID |
|---|---|---|
| chr4:70,344,431 A/G | 4:71,210,148 | rs34735123 |
| chr4:70,344,481 C/T | 4:71,210,198 | rs13101613 |
| chr4:70,344,604 A/T | 4:71,210,321 | rs67708234 |
| chr4:70,344,619 T/C | 4:71,210,336 | rs13134277 |
| chr4:70,345,253 A/G | 4:71,210,970 | rs13148087 |

These are region `4_71210148_71210970` of the reported set, taken from
`../../../4_SINGER_GenealogicalReconstruction_TSP/4_TSPRegions/3_GeneAssignment/variant_gene_assignments.pooled.tsv`.
A sixth TSP assigned to *MUC7*, rs28421149 (chr4:70,520,458), is a separate single-TSP region
(`4_71386175_71386175`) about 175 kb away and is not part of this locus.

## Scripts

- `eqtl_query.py` — eQTL Catalogue REST v2, for the two GTEx datasets of interest: minor salivary
  gland (QTD000276) and sun-exposed skin (QTD000316) → `assoc_rows.json`.
- `gtex_dyneqtl_v10.py` — GTEx portal API v2 `dyneqtl`, `datasetId=gtex_v10`, MUC7
  `ENSG00000171195.11`, for both tissues → `gtex_dyneqtl_v10.json`.
- `scan_tissues.py` — median TPM for MUC7 and its neighbours (CABS1, SMR3A, SMR3B, HTN1, HTN3)
  across all GTEx expression datasets, to establish where MUC7 is actually expressed.

All three query live APIs, so they need network access and their results depend on the database
version at the time of the call.

## Results

`results/` holds the tables as computed, from GTEx v10 throughout:

- `MUC7_eQTL_MSG_vs_skin.tsv` — the five TSPs in both tissues. Minor salivary gland p = 0.198
  (NES +0.044); sun-exposed skin p = 7.1e-8 to 1.3e-7 (NES -0.19).
- `MUC7_MSG_vs_skin_heterogeneity.tsv` — the two-tissue comparison per variant.
- `GTEx_v10_MUC7_MSG_significant_eQTLs.tsv`, `GTEx_v10_MUC7_median_expression.tsv`.

These come from live API calls, so they reflect the database at the time of the call.

## Software
Python 3 (standard library only). Network access to `www.ebi.ac.uk/eqtl/api/v2` and
`gtexportal.org/api/v2`.
