#!/usr/bin/env python3
"""Reduce the cobraa posterior grids to P(not yet coalesced by the human-chimp split).

The decode from ../3_Cobraa/ is a posterior over 64 coalescent-time intervals at every 1 kb
bin. This collapses it to a single column: the inverse cumulative posterior at the interval
boundary closest above 5.5 Mya, i.e. the probability that the two haplotypes of that individual
have NOT coalesced by the split.

MUTATION RATE. Mu used to convert the grid to years (N = theta / 4mu, years = T * 2N * gen) in
order to pick which column to take:

    mu = 1.27e-8 -> N = 19685 -> index 49 -> grid age 5,755,625 yr  (PO_5755625)   <- used here

Output: one BED per individual per chromosome, in panTro6 coordinates:
        CHROM  START  END  P(uncoalesced beyond the split)
        liftOver to hg38 happens next, in lift_clean_minimal_decode.sh.
"""
import os
import sys

import numpy as np
import pandas as pd

# ---- parameters (edit) ---------------------------------------------------------------------
# Decode output from ../3_Cobraa/, same parameter path
COBRAA_RESULTS = "../3_Cobraa/results"
PARAMS = ("D_64/b_100/spread1_0.075/spread2_100/muoverr_1.5/iterations30/thresh_1/"
          "thetafixed_0.001")
OUT_DIR = "minimal_decode/5.5Mya"

MU = 1.27e-8            # mutation rate per bp per generation
GEN = 24                # chimpanzee generation time, years
THETA_FIXED = 0.001     # as fitted in ../3_Cobraa/
SPLIT_AGE = 5.5e6       # human-chimp split, years

D, SPREAD_1, SPREAD_2 = 64, 0.075, 100      # must match the cobraa fit
EXPECTED_INDEX, EXPECTED_AGE = 49, 5755625  # guard: the grid point mu=1.27e-8 should select
# --------------------------------------------------------------------------------------------

SAMPLES = ["CEN_SAMEA4374772", "CEN_SAMEA4374778", "CEN_SAMEA4374785",
           "EAS_2003", "EAS_SAMEA4374790", "EAS_SAMEA4374797",
           "WES_1059", "WES_SAMEA2421542", "WES_SAMEA5204228"]
CHROMS = ["1", "2A", "2B", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14",
          "15", "16", "17", "18", "19", "20", "21", "22"]


def time_intervals(D, spread_1, spread_2):
    """cobraa's coalescent-time interval boundaries (scaled units)."""
    return np.array([0] + [spread_1 * np.exp((i / D) * np.log(1 + spread_2 / spread_1) - 1)
                           for i in range(D)])


os.makedirs(OUT_DIR, exist_ok=True)

T = time_intervals(D, SPREAD_1, SPREAD_2)
N = THETA_FIXED / (4 * MU)
T_years = T * 2 * N * GEN
idx = int(np.where(T_years > SPLIT_AGE)[0][0])
age_actual = int(T_years[idx])
assert (idx, age_actual) == (EXPECTED_INDEX, EXPECTED_AGE), (idx, age_actual)

colname = f"PO_{age_actual}"
print(f"mu={MU:.2e}  N={N:.1f}  cutoff index={idx}  grid age={age_actual}  (column {colname})")

n_done = 0
for popsam in SAMPLES:
    for chrom in CHROMS:
        infile = os.path.join(COBRAA_RESULTS, PARAMS, f"popsample_{popsam}",
                              f"chr{chrom}_decode.txt.gz")
        if not os.path.exists(infile):
            print(f"  MISSING {popsam} chr{chrom}", file=sys.stderr)
            continue
        arr = np.loadtxt(infile)
        position = arr[0, :]
        posterior = arr[1:, :]
        # inverse cumulative: P(coalescence time is at or beyond each interval)
        inv_cum = np.cumsum(posterior[::-1, :], axis=0)[::-1, :]
        pd.DataFrame({"CHROM": f"chr{chrom}",
                      "START": position.astype(int),
                      "END": (position + 1000).astype(int),
                      colname: inv_cum[idx, :]}).to_csv(
            f"{OUT_DIR}/{popsam}_chr{chrom}_minimal_decode.bed",
            sep="\t", index=False, header=False)
        n_done += 1
    print(f"  {popsam} done", flush=True)

print(f"wrote {n_done} beds to {OUT_DIR}")
