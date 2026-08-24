# 1_ChimpPSMC / 1_MappabilityMask — panTro6 mappability mask

Builds the mappability mask that `../2_Multihetsep/` intersects with each sample's callability
BED, so that regions where short reads cannot be placed uniquely are excluded from the
multihetsep.

Method is [SNPable](https://lh3lh3.users.sourceforge.net/snpable.shtml): extract every overlapping
k-mer of the reference (k = 150), realign them to the reference with `bwa aln -R 1000000 -O 3 -E 3`,
and score each position by the fraction of covering k-mers that map back uniquely and without
mismatch. That fraction is the stringency `r`.

## Files
- `run_snpable.sh` — the whole mask build, producing both stringencies used downstream:
  **r = 0.50** (main analysis) and **r = 0.75** (sensitivity analysis).
  ```
  ./run_snpable.sh <reference.fa> <outdir>
  ```
- `generate_bed.py` — converts a mask FASTA to BED, keeping positions scored 2 or 3.
  ```
  ./generate_bed.py PanTro6_final_mappability_mask.150.50.fa PanTro6_mappability_mask.150.50.bed
  ```
  
## Software
SNPable (`splitfa`, `gen_raw_mask.pl`, `gen_mask`); bwa; Python 3 with biopython.
