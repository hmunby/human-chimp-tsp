#!/usr/bin/env python3
"""Reorder a PLINK .fam so its rows follow the sample order of the BAM/CRAM list.

GATK's pedigree-aware calling requires the pedigree rows in the same order as the samples
in the merged VCF, which is the order of the BAM/CRAM list. Samples in the list that are
absent from the .fam are skipped.

    python3 reorder_fam.py chimp59_pedigree.fam bam_cram_chimp_samples.txt out.fam
"""
import sys

import pandas as pd

COLS = ["FAMILY_ID", "SAMPLE_ID", "PARENT1_ID", "PARENT2_ID", "SEX", "PHENOTYPE"]

fam_f, samples_f, out_f = sys.argv[1:4]

fam = pd.read_csv(fam_f, sep="\t", header=None, names=COLS)
samples = pd.read_csv(samples_f, sep="\t", header=None, names=["SAMPLE_ID", "FILE"])

order = [s for s in samples["SAMPLE_ID"] if s in set(fam["SAMPLE_ID"])]
out = fam.set_index("SAMPLE_ID").loc[order].reset_index()[COLS]
out.to_csv(out_f, sep="\t", header=False, index=False)
print(f"wrote {out_f}: {len(out)} of {len(fam)} pedigree rows, ordered by {samples_f}")
