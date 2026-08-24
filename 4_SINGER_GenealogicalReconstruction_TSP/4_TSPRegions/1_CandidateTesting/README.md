# 4_TSPRegions / 1_CandidateTesting — which candidates were testable, and why not

The TSP calling itself happens in `../../3_SINGER/`: `calc_tsp_all_shared.py` computes the
per-SNP support across a run's retained posterior ARGs, and `process_tsps.py` restricts to the
candidates and applies the support threshold. These scripts account for the denominator: of the
860 candidate SNPs, which ones reached the point of being testable at all.

## Result

| | candidates | |
|---|---|---|
| tested for TSP status | **802** | 93.3% |
| candidate lost in chimpanzee recalling | 27 | region excluded |
| candidate absent from the ARG | 31 | region analysed, position not in it |
| **total** | **860** | |

## Scripts

- `check_candidate_presence.py` — **operational pre-flight**, run after staging the merged region
  VCFs and before SINGER. A region whose candidate SNP did not survive recalling makes SINGER's
  downsample step fail and blocks aggregation, so those regions are excluded up front.
  `--write` appends them to the SINGER workdir's `misc_files/excluded_regions.txt`.

  It reports two things separately, because downsample matches on position alone while the
  correct test is chromosome **and** position. A region can be "rescued" by a position that
  coincides with a candidate on a different chromosome; those are flagged `masked-collision` and
  still excluded.

- `summarize_candidate_testing.py` — the table above. One row per candidate with the reason it
  was or was not tested. → `candidate_testing_summary.tsv`
- `summarize_lost_exclusions.py` — for each excluded region, why the candidate was lost:
  never called in chimpanzee, failed QC, or called and passing (so lost later, in merging or
  phasing). → `lost_exclusions_summary.tsv`

## Inputs to set

All three take `WORK`, the `../../3_SINGER/` working directory, at the top of the file.
`summarize_lost_exclusions.py` also takes `RECALL`, the `../../2_ChimpRecall/` working
directory, where it reads the per-chromosome `filtered/{chrom}/filtered.{chrom}_{pass,fail}.pos`
files.

The bundled `.tsv` files are the outputs from the paper run, kept as the reference result.

## Software
Python 3 (pandas); bcftools 1.20.
