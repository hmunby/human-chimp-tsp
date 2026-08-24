# Segmental-duplication track (chimp)

Generates the segmental-duplication mask used by the filtering Snakemake (`../Snakefile`, rule
`segdup`, stage 5). The goal is to remove segmental duplications that are **resolved (collapsed to
single copies) in the T2T chimpanzee assembly but not in the panTro6 reference** — those are a source
of false heterozygous / paralogous calls in panTro6-mapped data.

Method: [BISER](https://github.com/0xTCG/biser) compares two soft-masked assemblies and reports SDs both within and between them.

## Steps
| step | script | output |
|---|---|---|
| 0 (prep) | soft-mask both assemblies (see below) | masked FASTAs |
| 1 | `1_run_biser.sh` | `out/t2t_pantro6_seg_dups.txt` |
| 2 | `2_process_biser_out.py` | `processed_bed/…_resolved{,_no_chr}.bed` |
| 3 | `3_merge_sd_bed.sh` | `processed_bed/…_resolved{,_no_chr}.merged.bed` |

The final `t2t_pantro6_seg_dups_resolved_no_chr.merged.bed` (in `processed_bed/`) is the segdup
mask consumed by the filtering `segdup` rule (the Snakemake `dups` input in `../Snakefile`).

## Step 0 — masking (prerequisite)
- **panTro6**: UCSC soft-masked FASTA (download).
- **T2T (mPanTro3.hap1.cur.20231122)**: soft-masked by applying the assembly's own RepeatMasker
  annotation (`mPanTro3_v2.0.RepeatMasker_v1.1.hap1`, released with the T2T primate assembly) with
  `bedtools maskfasta`.

## Software
- BISER v1.4 (`pip install biser`)
- bedtools v2.29.1 (`maskfasta`, `merge`)
- Python 3 + pandas (`2_process_biser_out.py`)
