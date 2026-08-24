# 4_TSPRegions / 2_ClusteringLD — calling the TSP set and checking its regions

## The reported set: `aggregate_seeds.py`

The reported TSPs come from **three SINGER chains at Ne = 50,000 differing only in the MCMC
seed** (`../../3_SINGER`, runs `ne50000_seed1/2/3`), combined by pooling their posterior support.

    candidate SNP  ->  support rule  ->  single-linkage cluster at 10 kb  ->  regions

| rule | definition |
|---|---|
| **`pooled`** | **mean support over the three chains > 80** — reported |
| `strict` | passes independently in all three chains |
| `majority` | passes in at least two of three |
| `any` | passes in at least one |

`pooled` is the reported rule. Because each chain contributes exactly 100 trees, it is
identical to pooling all 300 posterior trees and asking whether more than 80 % show the
trans-species configuration, i.e. treating the three chains as the one longer chain they
are. `strict` is carried as the sensitivity set. Both are written by `aggregate_seeds.py`
and either can be passed to any downstream script.

Note the threshold convention: the count-based rules are sensitive to `>` versus `>=`
because a single chain can land on exactly 80, whereas a mean of three rarely does.
`--inclusive` switches it.

## Sensitivity runs: `cluster_tsps.py`

The same rule applied to a single chain, used for the Ne = 20,000 and Ne = 100,000 runs:

```bash
python cluster_tsps.py --run ../../3_SINGER --tag ne20k     # 111 TSPs, 49 regions
python cluster_tsps.py --run ../../3_SINGER_Ne100k --tag ne100k   # 133 TSPs, 70 regions
```

## Bundled TSP sets

Every downstream script defaults to `pooled` and takes a swappable parameter, so any of
these can be substituted without editing code:

| tag | TSPs | regions | what |
|---|---|---|---|
| `pooled` | 117 | 59 | **reported** — 3 chains at Ne = 50,000, pooled support |
| `strict` | 85 | 37 | same chains, must pass in all three |
| `ne20k` | 111 | 49 | single chain at Ne = 20,000 |
| `ne100k` | 133 | 70 | single chain at Ne = 100,000 |

`seed_support.tsv` records each candidate's support in all three chains plus every verdict.

## `within_region_ld.py`

Distance alone does not show that a region's TSPs sit on **one** trans-species haplotype.
This measures r² between every pair of TSPs in a region, in 1000 Genomes phase 3, within
AFR (n = 661), EUR (n = 503) and YRI (n = 108).

**Result for the reported set.** 12 of the 15 multi-TSP regions are a single tight block
(minimum pairwise r² ≥ 0.8 in every panel; median minimum 0.98). Three are not:

- **IGFBP7** (`4_57918296_57919705`) — two blocks, minimum r² 0.024 (AFR), 0.003 (EUR),
  0.038 (YRI) across the split against a mean of ~0.48. Two TSPs at 57,918,296–57,918,492
  and four at 57,919,221–57,919,705. Reported as its two blocks in Table 1.
- **MTRR** (`5_8017470_8024003`) — 0.858 in AFR, marginally below 0.8 in YRI only.
- **LINC01676** (`1_106125135_106131147`) — 0.808 in AFR, spanning 6 kb. This region is
  absent from `strict`.

## `region_genetic_span.py`

Physical span is not comparable across the genome; this converts each region's endpoints to
centimorgans on HapMap II (primary) and 1000G OMNI (sensitivity). Table 1.

No reported region exceeds **0.004 cM**, and the median multi-TSP region spans 0.0002 cM, so
every one is small enough in recombination units to be plausibly carried as one haplotype
since the human–chimpanzee split.

## Run

```bash
python aggregate_seeds.py                    # -> tsps.pooled.*, tsps.strict.*, ...
python within_region_ld.py                   # -> within_region_LD.{txt,tsv}
python region_genetic_span.py                # -> tsp_region_genetic_span*.tsv
# any of the above on another set, e.g.:
python within_region_ld.py --clusters tsp_clusters_10kb.strict.tsv --tsps tsps.strict.hg19.txt
```

## Not bundled

`resources/all_shared_snps_ann.txt` (shared-SNP annotation), the b37→hg38 chain, the
1000 Genomes phase 3 VCFs and panel, and the two genetic maps. See the stage README.

## Software

Python 3 (pandas, numpy); PLINK 1.9; bcftools; UCSC `liftOver`.
