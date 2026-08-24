# Allele ages (Relate, YRI)

The allele-age annotations are published Relate
estimates from Speidel et al. 2019 (*Nat. Genet.*), Zenodo record **3234689**
(https://zenodo.org/record/3234689), which provides per-population allele ages for 1000 Genomes
phase 3. We use **YRI** (from the AFR fileset). 

## Steps
1. Download the `allele_ages_YRI.RData` (YRI, in the AFR fileset) from Zenodo 3234689.
2. `reformat_relate_data.R <allele_ages_YRI.RData> <out.txt>` → a headerless tab file
   (`CHR, BP, BP, lower_age, upper_age, age_estimate`), then `bgzip` + `tabix -s1 -b2 -e3`.
   `age_estimate` is the **midpoint** of the lower/upper bounds (the point age estimate used downstream).
3. The driver adds INFO `LOWER_AGE_YRI/UPPER_AGE_YRI/AGE_ESTIMATE_YRI` (headers in `../misc_files/`).

Ages are in **generations**.

## Software
- R (for `reformat_relate_data.R`); bgzip / tabix.
