#!/usr/bin/env python3
"""Build the comprehensive per-SNP master table of the TSP set (117 SNPs for `pooled`).

This is a JOIN of pre-existing per-SNP tables -- nothing is recomputed here. Sources:

  1. THE TSP SET, already carrying most fields (coords in hg19+hg38, chimp panTro6 id,
     allele frequencies, CpG, BMAP, Roulette mutation rate, YRI age bounds, ARG support):
       ../../4_SINGER.../4_TSPRegions/2_ClusteringLD/tsps.pooled.hg19.txt

  2. GENE ASSIGNMENT + OPEN TARGETS evidence (join on hg38 coordinate):
       ../../4_SINGER.../4_TSPRegions/3_GeneAssignment/variant_gene_assignments.pooled.tsv
     -> rsID, nearest gene + distance, assigned gene + tier, VEP most-severe consequence,
        and the Open Targets molQTL / GWAS-L2G functional-evidence columns.

  3. phyloP (join on hg19 coordinate). PRIMARY column is the SAME source as the main figure:
       phylop/         = 36-eutherian-mammal EPO, HUMAN EXCLUDED
     SECONDARY comparison column:
       phylop_46way/   = UCSC 46-way placental, human INCLUDED
     Per-SNP scores read from shared_snps_phylop_chr{1..22}.txt (CHROM/POS in hg19).

Ages are stored in GENERATIONS in the source; years = generations x 29 (as MID_AGE_YRI_YEARS
in the source confirms). Midbranch = (lower + upper) / 2.

ARG support is the POOLED support: the mean of the three seed-replicate chains'
tsp_first_or_second_count. Because each chain contributes exactly 100 posterior genealogies,
that mean is the percentage of all 300 trees showing the trans-species configuration. It is
taken from `mean_support` and kept to one decimal -- `tsp_first_or_second_count` in the same
file is this value rounded to an integer, which would discard the fractional part that
distinguishes a pooled call from any single chain. For a single-chain set (ne20k, ne100k)
`mean_support` is that one chain's count, out of 100.

Column headers follow the manuscript convention (Table 1 / Supplementary Table): spelled out,
title case, units in parentheses.

Output: tsp_master_table.tsv  (one row per TSP; 117 for the pooled set, tab-separated)
"""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "../.."))
TR = f"{REPO}/4_SINGER_GenealogicalReconstruction_TSP/4_TSPRegions"

TSP = f"{TR}/2_ClusteringLD/tsps.pooled.hg19.txt"
GA = f"{TR}/3_GeneAssignment/variant_gene_assignments.pooled.tsv"
# 36-mammal EPO, human EXCLUDED -- the same track as the main phyloP figure. Large and not
# bundled; see 6_TSP_Analyses/snp_set_comparisons/README.md.
PHYLOP_MAIN = f"{REPO}/6_TSP_Analyses/snp_set_comparisons/phylop"
OUT = f"{HERE}/tsp_master_table.tsv"

# Optional overrides so the same table can be built for an alternative TSP set. Defaults are
# the reported pooled set; pass --tsps/--gene-assignments/--chains 1 for ne20k or ne100k.
import argparse as _argparse
_ap = _argparse.ArgumentParser()
_ap.add_argument("--tsps", default=TSP, help="clustered TSP table")
_ap.add_argument("--gene-assignments", default=None, help="variant_gene_assignments.tsv")
_ap.add_argument("--out", default=OUT)
_ap.add_argument("--chains", type=int, default=3,
                 help="posterior chains pooled in the supplied set (3 for pooled/strict, "
                      "1 for the single-chain ne20k/ne100k sets); sets the support header")
_a, _ = _ap.parse_known_args()
TSP = _a.tsps
OUT = _a.out
if _a.gene_assignments:
    GA = _a.gene_assignments

GEN_TO_YR = 29  # generation time (years) used throughout the manuscript

# Roulette reports a scaled rate, not a per-generation one. This is the same factor Figure 1
# applies to put its mutation-rate panel in per-site-per-generation units, so the table and
# the figure are on one scale. Roulette's own mean over 1000GP MAF>=5% sites is 0.354944 on
# the raw scale, i.e. 3.60e-8 per site per generation.
ROULETTE_SCALE = 1.015e-7

# Support header. The pooled set averages three 100-tree chains (300 genealogies); the
# single-chain sensitivity sets have 100. Set by --chains so the header never misstates it.
SUPPORT_COL = ("ARG support for TSP (% of 300 posterior genealogies, 3 chains x 100)"
               if _a.chains == 3 else
               f"ARG support for TSP (% of {_a.chains * 100} posterior genealogies)")

# --------------------------------------------------------------------------------- load the TSP set
t = pd.read_csv(TSP, sep="\t", low_memory=False)
# row count is whatever the supplied TSP set holds (117 for the reported pooled set)
print(f"TSP set: {len(t)} rows from {TSP}")

# chimp panTro6 coordinate is encoded in CHIMP_ID = "chrom_pos_ref_alt"
chimp = t["CHIMP_ID"].str.split("_", expand=True)
t["chimp_chrom"] = chimp[0]
t["chimp_pos"] = pd.to_numeric(chimp[1], errors="coerce").astype("Int64")

# ancestral allele: AA field is "A|||" style; take the first token. Lowercase in the 1000G
# ancestral call denotes a low-confidence assignment -> flag it, then upper-case the base.
aa_raw = t["AA"].fillna("").str.split("|").str[0]
t["anc_lowconf"] = aa_raw.str.islower()
t["anc"] = aa_raw.str.upper().replace("", np.nan)

# numeric coercions for age fields ("." -> NaN where Relate has no age)
for c in ["LOWER_AGE_YRI", "UPPER_AGE_YRI", "MID_AGE_YRI_YEARS", "BMAP", "MR", "AF", "CHIMP_AF"]:
    t[c] = pd.to_numeric(t[c], errors="coerce")

t["mid_age_gen"] = (t["LOWER_AGE_YRI"] + t["UPPER_AGE_YRI"]) / 2.0
t["lower_age_yr"] = t["LOWER_AGE_YRI"] * GEN_TO_YR
t["upper_age_yr"] = t["UPPER_AGE_YRI"] * GEN_TO_YR
t["mid_age_yr"] = t["mid_age_gen"] * GEN_TO_YR
# sanity: source MID_AGE_YRI_YEARS should equal our mid_age_yr where both present
chk = t.dropna(subset=["MID_AGE_YRI_YEARS", "mid_age_yr"])
assert np.allclose(chk["MID_AGE_YRI_YEARS"], chk["mid_age_yr"], rtol=1e-4), "midbranch-year mismatch"

t["cpg_label"] = t["CPG"].map({1: "CpG", 0: "non-CpG"}).fillna("NA")
t["hg38_key"] = t["CHROM_HG38"].astype(str) + "_" + t["POS_HG38"].astype(str)

# ----------------------------------------------------------------------------- gene assignment + OT
v = pd.read_csv(GA, sep="\t", dtype=str, low_memory=False)
v["hg38_key"] = v["CHROM_HG38"].astype(str) + "_" + v["POSITION_HG38"].astype(str)
vcols = ["hg38_key", "rsid", "most_severe", "nearest_gene", "nearest_distance",
         "assigned_gene", "assignment_tier", "in_ot_index", "n_credible_sets",
         "molqtl_details", "gwas_l2g_genes"]
v = v[vcols].drop_duplicates("hg38_key")
m = t.merge(v, on="hg38_key", how="left", validate="one_to_one")
assert m["rsid"].notna().sum() >= 1, "gene-assignment join produced no rsIDs"
n_unmatched = m["assignment_tier"].isna().sum()
print(f"gene-assignment join: {len(m)-n_unmatched}/{len(m)} matched on hg38 coord")

# molQTL type(s) and tissue(s): unique, order-preserving, across ALL colocalising credible sets.
# molqtl_details is ";"-separated entries of the form "gene|type|tissue|pip=X|method".
def _uniq_field(details, idx):
    if pd.isna(details):
        return np.nan
    vals = []
    for entry in str(details).split(";"):
        parts = entry.split("|")
        if len(parts) > idx:
            val = parts[idx].strip()
            if val and val not in vals:
                vals.append(val)
    # "; " separator (tissue names themselves can contain commas)
    return "; ".join(vals) if vals else np.nan

m["molqtl_types"] = m["molqtl_details"].apply(lambda x: _uniq_field(x, 1))
m["molqtl_tissues"] = m["molqtl_details"].apply(lambda x: _uniq_field(x, 2))

# ----------------------------------------------------------------------------- phyloP (hg19 join)
def load_phylop(dirpath):
    parts = [pd.read_csv(f"{dirpath}/shared_snps_phylop_chr{c}.txt", sep="\t") for c in range(1, 23)]
    p = pd.concat(parts, ignore_index=True)
    p["hg19_key"] = p["CHROM"].astype(str) + "_" + p["POS"].astype(str)
    return p.drop_duplicates("hg19_key").set_index("hg19_key")["PHYLOP_SCORE"]

m["hg19_key"] = m["CHROM"].astype(str) + "_" + m["POS"].astype(str)
ph_main = load_phylop(PHYLOP_MAIN)
m["phylop_main"] = m["hg19_key"].map(ph_main)
print(f"phyloP (36-mammal EPO, human-excluded, main): {m['phylop_main'].notna().sum()}/{len(m)} matched")

# ----------------------------------------------------------------------------- assemble output
out = pd.DataFrame({
    "rsID": m["rsid"],
    "Chromosome (hg19)": m["CHROM"],
    "Position (hg19)": m["POS"],
    "Chromosome (hg38)": m["CHROM_HG38"],
    "Position (hg38)": m["POS_HG38"],
    "Chimp chromosome (panTro6)": m["chimp_chrom"],
    "Chimp position (panTro6)": m["chimp_pos"],
    "Reference allele": m["REF"],
    "Alternate allele": m["ALT"],
    "Ancestral allele": m["anc"],
    "Allele frequency (1000 Genomes)": m["AF"].round(4),
    "Allele frequency (chimpanzee)": m["CHIMP_AF"].round(4),
    "CpG status": m["cpg_label"],
    "Background-selection B value (BMAP)": m["BMAP"].astype("Int64"),
    "Mutation rate, Roulette (per site per generation)":
        (pd.to_numeric(m["MR"], errors="coerce") * ROULETTE_SCALE).map(
            lambda x: f"{x:.3e}" if pd.notna(x) else ""),
    "phyloP (36-mammal EPO, human excluded)": m["phylop_main"].round(3).fillna("NA"),
    "Variant annotation (snpEff)": m["Annotation"],
    "Most severe consequence (VEP)": m["most_severe"],
    "Nearest gene": m["nearest_gene"],
    "Distance to nearest gene (bp)": pd.to_numeric(m["nearest_distance"], errors="coerce").astype("Int64"),
    "Assigned gene": m["assigned_gene"],
    "Assignment tier": m["assignment_tier"],
    "In Open Targets index": m["in_ot_index"],
    "Number of credible sets": pd.to_numeric(m["n_credible_sets"], errors="coerce").astype("Int64"),
    "molQTL type(s)": m["molqtl_types"],
    "molQTL tissue(s)": m["molqtl_tissues"],
    "GWAS L2G effector genes": m["gwas_l2g_genes"],
    "Lower age bound, YRI (years)": m["lower_age_yr"].round(0).astype("Int64"),
    "Upper age bound, YRI (years)": m["upper_age_yr"].round(0).astype("Int64"),
    "Midbranch age, YRI (years)": m["mid_age_yr"].round(0).astype("Int64"),
    SUPPORT_COL: pd.to_numeric(m["mean_support"], errors="coerce").round(1),
})

out = out.sort_values(["Chromosome (hg19)", "Position (hg19)"],
                      key=lambda s: pd.to_numeric(s, errors="coerce") if s.name == "Chromosome (hg19)" else s
                      ).reset_index(drop=True)
out.to_csv(OUT, sep="\t", index=False)
print(f"\nwrote {OUT}  ({len(out)} rows x {out.shape[1]} columns)")
print("columns:")
for c in out.columns:
    print("  -", c)
