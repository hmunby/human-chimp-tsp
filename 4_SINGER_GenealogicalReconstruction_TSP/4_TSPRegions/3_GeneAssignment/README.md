# 4_TSPRegions / 3_GeneAssignment — one gene per TSP region

`opentargets_assign.py` assigns each of the 59 TSP regions a single gene, preferring direct
functional variant-to-gene evidence over proximity.

## The rule

Applied per region, not per SNP, so that all TSPs in a region inherit one decision and the
per-variant and per-region tables can never disagree:

| tier | evidence | regions |
|---|---|---|
| `l2g` | GWAS locus-to-gene effector prediction, score ≥ 0.5 (L2G already integrates molQTL colocalisation, distance, chromatin and VEP) | 14 |
| `molqtl` | membership of a molecular-QTL 95% credible set (eQTL, sQTL, pQTL, tuQTL, edQTL, sceQTL, scsQTL). Where a region has several molQTL genes, the winner is chosen on breadth: tissues, then modalities, then credible sets | 14 |
| `nearest` | `bedtools closest` against GENCODE v49 protein-coding autosomal genes | 31 |

59 regions → **55 unique genes**.

Assigning per region on aggregated evidence matters: a single near-zero-PIP eQTL can outrank the
region's real signal if each SNP is decided alone: deciding per SNP sends the GYPB region to the
*GUSBP5* pseudogene and MTRR to an unnamed ENSG model.

## Coverage

Of the 117 TSPs, 75 are in the OpenTargets variant index, 53 have a colocalising molQTL and 32 a
GWAS L2G link. Variants absent from the index (roughly 40% of arbitrary SNPs) legitimately have no
functional evidence and fall through to the nearest-gene tier. Every variant's PIP within each
credible set and that set's confidence are recorded, so a PIP threshold can be applied afterwards
without re-querying.

## Run

```
python opentargets_assign.py
```

Queries the OpenTargets Platform GraphQL API in batches of 20 variants and caches every response
in `ot_variant_cache.json`, so re-runs are offline. Variant ids are `CHR_POS_REF_ALT` on
**GRCh38**, taken from the `CHROM_HG38`/`POSITION_HG38` columns.

Outputs, one pair per TSP-set tag:
- `region_gene_assignments.pooled.tsv` — 59 rows, the foreground gene set.
- `variant_gene_assignments.pooled.tsv` — 117 rows, the per-variant evidence behind them.

The `.strict`, `.ne20k` and `.ne100k` pairs are the sensitivity sets. All are bundled as the
paper's result. The gene-set enrichment built on the foreground set, and the ascertainment-matched
null loci it needs, live in `../../../6_TSP_Analyses/enrichment/geneset/`.

## Inputs to set
- `SNP_TABLE` — `../2_ClusteringLD/tsps.pooled.for_gene_assignment.txt` (pre-filled).
- `GENCODE_BED` — GENCODE v49 protein-coding autosomal genes as BED with gene names; not bundled.

## Note on the API
OpenTargets Genetics has been merged into the Platform, so the endpoint is
`api.platform.opentargets.org/api/v4/graphql`. The schema this script targets was current in
June 2026; the cached JSON is what the paper uses. 

## Software
Python 3 (pandas); bedtools; network access for the OpenTargets API.
