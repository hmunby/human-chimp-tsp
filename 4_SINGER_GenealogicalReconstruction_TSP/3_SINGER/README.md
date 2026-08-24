# 4 / 3_SINGER — SINGER genealogical reconstruction and the per-run TSP calls

Runs SINGER per recall region on the merged human+chimp VCFs from `../2_ChimpRecall/`, converts
the ARGs to tskit, computes per-SNP TSP support, and calls the TSPs for each run.

One Snakefile covers every configuration in the paper:

| run | Ne | seed | |
|---|---|---|---|
| `ne50000_seed1` | 50,000 | 1 | the three chains of the **reported** set |
| `ne50000_seed2` | 50,000 | 2 | |
| `ne50000_seed3` | 50,000 | 3 | |
| `ne20000_seed1` | 20,000 | 1 | Ne sensitivity |
| `ne100000_seed1` | 100,000 | 1 | Ne sensitivity |

They differ only in `-Ne` and `-seed`. Everything upstream is shared: the merged region VCFs and
the downsampled 50-human panel are built once and consumed byte-identically by all five, so no run
can perturb another.

The reported set is not produced here. Each run yields its own single-chain call; the three
Ne = 50,000 chains are pooled by `../4_TSPRegions/2_ClusteringLD/aggregate_seeds.py`.

## Setup

Copy the per-region merged VCFs from `../2_ChimpRecall/` into `merged_regions/` and index them
(see the Snakefile header). Regions are the intersection of what is staged with the 250 kb target
list, minus `misc_files/excluded_regions.txt` — regions whose candidate SNP was lost in chimp
recalling, so SINGER has nothing to anchor on.

## Run

```bash
snakemake --profile .                                             # all five
snakemake --profile . --config runs=ne20000_seed1,ne50000_seed1   # a subset
snakemake --profile . --config archive_arg=1                      # tar each region's raw ARG
```

## Flow (rules, in order)

- **downsample** — subset the human samples to 50 per region, sampling from each diploid genotype
  class at the candidate SNP in the proportions of the full 1000GP panel, so the focal allele
  frequency is preserved (`scripts/downsample_human_samples.py`). Run once, shared by all runs.
- **run_singer** — `singer_master -Ne <ne> -seed <seed> -m 1.25e-8 -n 167 -thin 60` per region.
- **convert_to_tskit** — ARG samples → tskit `.trees` for the retained posterior samples 67–166
  (0–66 are burn-in). Uses a patched converter that sorts the mutation table and computes mutation
  parents; SINGER's bundled one omits both and fails on strict tskit.
- **all_shared_metrics** — per-SNP support, scored at **every shared
  SNP** in the region rather than only the candidates, so the support landscape around each locus
  is available if needed. 
- **archive_arg** — optional, tars a region's raw ARG once its metrics exist.
- **process_tsps** — that run's TSP call:

      candidate SNP  ->  support > 80%  ->  TSPs

  Restricted to the original 860 candidates, Outputs `results/<run>/candidate_support.txt` (every candidate with its
  support, passed or not — the input to `aggregate_seeds.py`), `results/<run>/tsps.txt` and
  `results/<run>/tsps.hg38.bed`.

## Inputs to set (top of `Snakefile`, `resources/`)

- `SHARED_SITES` — shared-SNP sites from stage 3 (relative path pre-filled).
- `CANDIDATES` — the 860 candidate SNPs (bundled in `misc_files/`).
- `SNP_INFO` — annotated shared-SNP table (`all_shared_snps_ann.txt`), not bundled.
- `CHAIN` — b37 → hg38 UCSC liftOver chain, not bundled.

Resource requests scale with each region's VCF size and with Ne; SINGER holds the ARG in memory,
and Ne = 100,000 needs roughly 3× the resident set of Ne = 20,000. Dropping a
`misc_files/region_hours.json` of measured per-region wall clock lets the time request be sized
per region instead of every job asking the same ceiling.

## A note on file count

Five runs × ~433 regions × ~500 small files is ~1M files: SINGER writes four text files per
posterior sample and the converter one `.trees` each. That is fine on a filesystem sized for it.
If yours objects, run the configurations one at a time or set `archive_arg=1`.

## Software

SINGER (`singer_master`); tskit; Python 3 (numpy, pandas); bcftools 1.20; UCSC liftOver;
Snakemake 7.
