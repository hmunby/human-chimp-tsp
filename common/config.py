"""Central configuration for the paper reproducibility scripts.

Input-data paths live here. Inputs that are too large to bundle are expected under `resources/` beside this file; set PAPER_DATA_ROOT to put
them elsewhere. The small derived tables the scripts actually read are bundled.

    export PAPER_DATA_ROOT=/path/to/resources   # optional override
"""
import os

# --- roots -----------------------------------------------------------------------------------
PROJECT_ROOT = os.environ.get(
    "PAPER_DATA_ROOT",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources"),
)
ANALYSES = PROJECT_ROOT

# output directory for regenerated figures (override with PAPER_FIG_OUTDIR)
FIG_OUTDIR = os.environ.get("PAPER_FIG_OUTDIR", os.path.join(os.path.dirname(__file__), "figures"))

# --- shared constants ------------------------------------------------------------------------
GEN_YEARS = 29            # years per generation (human)
ROULETTE_SCALE = 1.015e-7  # Roulette relative MR * ROULETTE_SCALE = per-site-per-generation rate
TGP_MEAN_MR = 0.354944     # mean Roulette MR of 1000GP MAF>=5% SNPs (Panel B reference line)
HC_SPLIT_MYA = (5, 6)      # human-chimp split band shaded in the age-distribution panels

# --- Figure 1: age distribution + mutation rate ----------------------------------------------
# Panel A curves are precomputed for each set (MHC-included and MHC-excluded). Panel B reads the
# shared-SNP annotation table and applies the MHC filter itself.
_AGE_OUT = os.path.join(ANALYSES, "age_snp_uncoal_analysis/outputs")

MHC_B37 = ("6", 28_477_797, 33_448_354)   # MHC in b37/hg19, as used throughout

FIG1_TGP_EXPECTED = os.path.join(_AGE_OUT, "tgp_hg38_expected_snps.step_size_50.txt")
FIG1_SHARED_EXPECTED = os.path.join(_AGE_OUT, "shared_expected_snps.step_size_50.txt")

FIG1_TGP_EXPECTED_NOMHC = os.path.join(_AGE_OUT, "tgp_hg38_no_mhc_expected_snps.step_size_50.txt")
FIG1_SHARED_EXPECTED_NOMHC = os.path.join(_AGE_OUT, "shared_no_mhc_expected_snps.step_size_50.txt")
# Per-SNP shared-SNP annotation (hg19): CHROM, POS, MR, LOWER_AGE_YRI, UPPER_AGE_YRI, ...
FIG1_SHARED_ANN = os.path.join(PROJECT_ROOT, "all_shared_snps_ann.txt")
