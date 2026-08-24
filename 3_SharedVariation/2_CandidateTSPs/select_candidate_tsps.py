#!/usr/bin/env python3
"""Select the candidate trans-species polymorphisms from the shared-variant set.

A shared SNP is a candidate if its midpoint YRI Relate age is >= 4 Mya, excluding the MHC
(b37 chr6:28,477,797-33,448,354). Midpoint age is the mean of the LOWER_AGE_YRI and
UPPER_AGE_YRI bounds, in generations, times 29. This yields the 860 candidates that stage 4
recalls in chimpanzee and reconstructs with SINGER.

    python select_candidate_tsps.py

Input:  the shared-variant table from ../1_Intersect/ (bcftools query output, no header).
Output: candidate_tsps.chrom_pos.txt  (CHROM<TAB>POS).

Python 3 + pandas.
"""
import pandas as pd

# ---- parameters (edit) --------------------------------------------------------------------------
SHARED_TABLE = "../1_Intersect/paralog_filtered/shared_variants.SD_CNV_DUP_filtered.table.txt"
OUT = "candidate_tsps.chrom_pos.txt"
GEN_YEARS = 29                        # years per generation
MIN_AGE_YEARS = 4_000_000             # candidate threshold: midpoint YRI age >= 4 Mya
MHC = ("6", 28_477_797, 33_448_354)   # MHC region excluded (b37 / hg19 chr6)
# ------------------------------------------------------------------------------------------------

# columns of the intersect variant table (see ../1_Intersect/Snakefile rule output_variant_table)
COLS = ["CHROM", "POS", "REF", "ALT", "AF", "CHIMP_ID", "CHIMP_AF", "AA", "CPG",
        "LOWER_AGE_YRI", "UPPER_AGE_YRI", "ANN"]
df = pd.read_csv(SHARED_TABLE, sep="\t", header=None, names=COLS)

# midpoint YRI age, generations -> years
mid_age_years = df[["LOWER_AGE_YRI", "UPPER_AGE_YRI"]].mean(axis=1) * GEN_YEARS

# exclude the MHC region
mhc_chrom, mhc_start, mhc_end = MHC
in_mhc = (df["CHROM"].astype(str) == mhc_chrom) & df["POS"].between(mhc_start, mhc_end)

candidates = df[(mid_age_years >= MIN_AGE_YEARS) & ~in_mhc].sort_values(["CHROM", "POS"])
candidates[["CHROM", "POS"]].to_csv(OUT, sep="\t", header=False, index=False)
print(f"candidate TSPs (midpoint YRI age >= {MIN_AGE_YEARS/1e6:.0f} Mya, MHC excluded): {len(candidates)}")
