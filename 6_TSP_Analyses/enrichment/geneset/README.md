# enrichment / geneset — GO, KEGG, Reactome and Hallmark, against an ascertainment-matched null

The TSP loci are not a random sample of the genome: they are ascertained in high-BMAP (weak
background selection) regions, and the gene-assignment rule (L2G → supported molQTL → nearest)
has its own structural biases toward large, well-studied genes. A whole-genome background, which
is what a default DAVID or Enrichr run gives, therefore reports categories that reflect the
ascertainment rather than the biology. `glycoprotein` is the worked example: strongly enriched
against a naive background, not enriched at all once ascertainment is matched.

Every category is therefore tested **twice**, and the pair is the result:

| test | background | reads as |
|---|---|---|
| `unmatched` | all autosomal protein-coding genes, hypergeometric | what a default whole-genome background gives |
| `matched` | 2,000 BMAP + span + SNP-count matched null gene sets, empirical | what survives ascertainment |

A category significant against the unmatched but **not** the matched null is an ascertainment
artefact. One significant against **both** is a candidate real signal.

## Components

- `build_background.py` — samples null loci matched to the foreground on BMAP (± 20),
  SNP-cluster span and SNP count, then pushes them through the *identical* `assign_one_gene`
  rule imported from `../../../4_SINGER_GenealogicalReconstruction_TSP/4_TSPRegions/3_GeneAssignment/opentargets_assign.py`.
  Two stages: stage 1 samples and validates the matching offline; stage 2 queries OpenTargets
  (cached) and emits the background universe and per-replicate null gene sets. Review stage 1
  before running stage 2 — the API stage is the expensive one.
- `enrichment_scan.py` — the two-way scan. `--regions` selects the foreground, `--null` the
  null gene sets, `--min-k` the minimum foreground genes for a category to be tested.
- `reassign_null.py` — re-runs the assignment rule over an existing null manifest, for when the
  rule changes but the sampled loci do not.
- `glycoprotein_enrichment_test.py`, `glycoprotein_enrichment_test.memglyco.py` — the
  single-category permutation tests for the UniProt glycoprotein (KW-0325) and membrane
  glycoprotein (KW-0325 ∩ KW-0472) sets. The distribution figures for these live in
  `../glycoprotein/`.

## Run

```
python enrichment_scan.py --min-k 3 --out enrichment_scan_results.MINK3.tsv
python glycoprotein_enrichment_test.py
```

To rebuild the null from scratch (needs network for the OpenTargets stage):

```
python build_background.py --stage 1     # sample and validate matching, offline
python build_background.py --stage 2     # assign genes, emits null_gene_sets.tsv
```

## Bundled

| file | what |
|---|---|
| `genesets/*.gmt` | GO BP/CC/MF 2023, KEGG 2021 Human, Reactome 2022, MSigDB Hallmark 2020 |
| `null_gene_sets.tsv` | 2,000 matched null gene sets (the published null) |
| `background_genes.txt` | the null gene universe, 12,793 genes |
| `enrichment_scan_results.tsv` | the k ≥ 2 scan |
| `enrichment_scan_results.MINK3.tsv` | the k ≥ 3 scan, the reported one |
| `glycoprotein_enrichment_results.tsv`, `membrane_glycoprotein_enrichment_results.tsv` | the two single-category tests |

## Not bundled

`resources/gencode.v49.annotation.genes.protein_coding.autosomes.gene_names.bed` (bed file derived from restricting GFF to protein coding genes), and `resources/null_locus_assignments.tsv` / `null_loci_manifest.tsv` (7 MB and 10 MB;
the per-locus null sampling record, needed only to resample or re-assign).

## Software

Python 3 (pandas, numpy, scipy); bedtools; network access only for `build_background.py` stage 2.
