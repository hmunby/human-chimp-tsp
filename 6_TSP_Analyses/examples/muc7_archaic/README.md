# MUC7 TSPs versus the archaic haplotype

*MUC7* carries a deeply divergent haplotype in African populations that Xu et al. (2017,
*Mol Biol Evol* 34:2704–2715) attribute to introgression from an unidentified archaic hominin.
The results text states that the TSPs we identify lie close to, but are distinct from and not in
linkage disequilibrium with, that haplotype. `muc7_archaic_ld.py` measures it.

## Result

The archaic haplotype is Xu et al.'s 44 group-E tag SNPs, spanning chr4:71,337,964–71,348,157
(hg19). The MUC7 TSPs flank it without overlapping: the five-TSP region at 71,210,148–71,210,970
sits ~127 kb proximal, the single-TSP region at 71,386,175 ~38 kb distal.

Maximum r² between each TSP and any of the 44 tag SNPs, over all 264 pairs:

| TSP (hg19) | AFR (n = 661) | all 1000GP |
|---|---|---|
| chr4:71,210,148 | 0.021 | 0.073 |
| chr4:71,210,198 | 0.021 | 0.073 |
| chr4:71,210,321 | 0.021 | 0.073 |
| chr4:71,210,336 | 0.021 | 0.073 |
| chr4:71,210,970 | 0.021 | 0.073 |
| chr4:71,386,175 | 0.025 | 0.073 |

Median r² is 0.005 in AFR and 0.003 overall. 

## Run

```bash
python muc7_archaic_ld.py        # -> muc7_archaic_ld.tsv, one row per (TSP, tag SNP) pair
```

## Files

- `xu_2017_data/haplotype_snps.txt` — the 44 group-E tag SNPs, bundled.
- `muc7_archaic_ld.tsv` — every (TSP, tag SNP) pair with its r², in both panels.

## Software

Python 3 (pandas); bcftools; PLINK 1.9.
