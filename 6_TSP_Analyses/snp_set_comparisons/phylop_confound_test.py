#!/usr/bin/env python3
"""Is the shared-SNP phyloP deficit driven by the chimpanzee reference-allele artefact?

Chimpanzee is present in the 36-eutherian EPO alignment (human is masked to N), so a site that
is variable in chimpanzee can have its phyloP depressed whenever the chimpanzee reference
happens to carry the derived allele. That probability rises with the chimpanzee derived-allele
frequency, so the artefact predicts phyloP should DECLINE as chimp DAF rises.

Competing explanation: shared SNPs are enriched for recurrent mutation at hypermutable (CpG)
sites, which evolve fast across mammals and score low regardless of the artefact.

Stratifying shared SNPs by chimp DAF and by CpG status separates the two, and the TSPs are
placed on both axes for comparison.

Inputs:  phylop/ per-chromosome score tables from ./Snakefile; the shared-SNP annotation table
         for CHIMP_AF / CpG / ancestral allele; misc_files/ TSP list.
Output:  printed to stdout.
"""
import os

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))

# ---- inputs (edit) -------------------------------------------------------------------------
# Shared-SNP annotation table (hg19): CHROM, POS, REF, ALT, AA, AF, CHIMP_AF, CPG. Same file
# the TSP calling uses, and hg19 like the phyloP tables. Large, not bundled -- see README.
SHARED_TBL = os.path.join(HERE, "resources/all_shared_snps_ann.txt")
PHYLOP_DIR = os.path.join(HERE, "phylop")
TSP_FILE = os.path.join(HERE, "misc_files/tsp_snps_pooled.hg19.txt")
# --------------------------------------------------------------------------------------------

ph = pd.concat([pd.read_csv(f"{PHYLOP_DIR}/shared_snps_phylop_chr{c}.txt", sep="\t")
                for c in range(1, 23)])
ph = ph.dropna(subset=["PHYLOP_SCORE"])

cols = ["CHROM", "POS", "CHIMP_AF", "AF", "CPG", "AA", "REF", "ALT"]
tbl = pd.read_csv(SHARED_TBL, sep="\t", usecols=cols, low_memory=False)
# Ensembl writes the ancestral allele as "T|||", lowercase where the call is low-confidence.
tbl["ANC"] = tbl["AA"].astype(str).str[0].str.upper()
tbl = tbl[tbl["ANC"].isin(list("ACGT"))].copy()
for c in ("CHROM", "POS"):
    tbl[c] = pd.to_numeric(tbl[c], errors="coerce")
tbl["CHIMP_AF"] = pd.to_numeric(tbl["CHIMP_AF"], errors="coerce")
tbl = tbl.dropna(subset=["CHROM", "POS", "CHIMP_AF"]).astype({"CHROM": int, "POS": int})

d = pd.merge(tbl, ph, on=["CHROM", "POS"], how="inner")
print(f"shared SNPs with phyloP and CHIMP_AF: {len(d):,}")

# chimpanzee DERIVED allele frequency: CHIMP_AF is the ALT frequency, so flip if ALT is ancestral
d["chimp_daf"] = np.where(d["ALT"] == d["ANC"], 1 - d["CHIMP_AF"], d["CHIMP_AF"])
print(f"  chimp DAF: median {d['chimp_daf'].median():.3f}  "
      f"IQR [{d['chimp_daf'].quantile(.25):.3f}, {d['chimp_daf'].quantile(.75):.3f}]")

print("\n=== phyloP vs chimpanzee derived-allele frequency (artefact predicts a DECLINE) ===")
d["afbin"] = pd.cut(d["chimp_daf"], [0, .05, .10, .20, .30, .50, .70, 1.01], right=False)
print(f"  {'chimp DAF':>14} {'n':>8} {'median phyloP':>14} {'mean':>8}")
for b, g in d.groupby("afbin", observed=True):
    print(f"  {str(b):>14} {len(g):>8,} {g['PHYLOP_SCORE'].median():>+14.3f} "
          f"{g['PHYLOP_SCORE'].mean():>+8.3f}")
rho, p = spearmanr(d["chimp_daf"], d["PHYLOP_SCORE"])
print(f"  Spearman rho = {rho:+.4f}  p = {p:.2e}")

print("\n=== phyloP by CpG status (hypermutability explanation) ===")
d["CPG"] = pd.to_numeric(d["CPG"], errors="coerce")
print(f"  {'CpG':>6} {'n':>8} {'median phyloP':>14} {'mean':>8}")
for v, g in d.groupby("CPG", observed=True):
    print(f"  {int(v):>6} {len(g):>8,} {g['PHYLOP_SCORE'].median():>+14.3f} "
          f"{g['PHYLOP_SCORE'].mean():>+8.3f}")
cpg, noncpg = d[d.CPG == 1]["PHYLOP_SCORE"], d[d.CPG == 0]["PHYLOP_SCORE"]
print(f"  CpG vs non-CpG: MWU p = {mannwhitneyu(cpg, noncpg, alternative='two-sided')[1]:.2e}  "
      f"Δmedian = {cpg.median() - noncpg.median():+.3f}")

print("\n=== joint: median phyloP by CpG x chimp DAF ===")
print(d.pivot_table(index="afbin", columns="CPG", values="PHYLOP_SCORE",
                    aggfunc="median", observed=True).to_string())

# where do the TSPs sit on these axes?
tsp = pd.read_csv(TSP_FILE, sep="\t")
tsp.columns = [c.upper() for c in tsp.columns]
t = pd.merge(tsp[["CHROM", "POS"]], d, on=["CHROM", "POS"], how="inner")
print(f"\n=== TSPs (n={len(t)}) ===")
print(f"  chimp DAF median {t['chimp_daf'].median():.3f}   CpG fraction {t['CPG'].mean():.3f}"
      f"   median phyloP {t['PHYLOP_SCORE'].median():+.3f}")
print(f"  all shared: chimp DAF median {d['chimp_daf'].median():.3f}   "
      f"CpG fraction {d['CPG'].mean():.3f}   median phyloP {d['PHYLOP_SCORE'].median():+.3f}")
