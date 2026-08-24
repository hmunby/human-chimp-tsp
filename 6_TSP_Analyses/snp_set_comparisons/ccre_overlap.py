#!/usr/bin/env python3
"""Overlap of each SNP set with ENCODE SCREEN candidate cis-regulatory elements.

The same four nested sets used everywhere else in this directory (1000GP MAF >= 5% >
all shared > old shared >= 4 Mya > TSP), scored as the percentage of SNPs falling in any
cCRE. Positions and the registry are both GRCh38.

The reported registry is V4 (2,348,854 elements, 21.2% of the autosomes): TSP 15.4%,
old shared 18.5%, all shared 20.7%, 1000GP 22.3%. V3 is selectable.

Usage:
    python ccre_overlap.py                      # reported: pooled TSPs, V4 registry
    python ccre_overlap.py --registry v3         # against the 2020 registry
    python ccre_overlap.py --tsp-bed misc_files/tsp_snps_strict.hg38.bed
"""
import argparse
import os
import subprocess

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))

# ---- inputs (edit) -------------------------------------------------------------------------
# ENCODE SCREEN Registry of cCREs, GRCh38. Not bundled (63-129 MB) -- see README.
#   v4  Moore et al. Nature (2026), doi:10.1038/s41586-025-09909-9 -- the reported registry
#   v3  Moore et al. Nature 583:699-710 (2020)
REGISTRY = {
    "v4": "resources/functional_annotation_datasets"
          "ENCODE_SCREEN_CCRES/2026/GRCh38-cCREs.bed",
    "v3": "resources/functional_annotation_datasets"
          "ENCODE_SCREEN_CCRES/GRCh38-cCREs.bed",
}

# hg38 SNP-position BEDs. The two large ones are not bundled -- see README.
TGP_BED = os.path.join(HERE, "resources/tgp_snps_hg38.bed")                          # 153 MB
SHARED_BED = os.path.join(HERE, "resources/shared_snps_chrom_pos_hg38.sorted.bed")   # 1.6 MB
OLD_SHARED_BED = os.path.join(HERE, "misc_files/old_shared_snps_860.hg38.sorted.bed")
TSP_BED = os.path.join(HERE, "misc_files/tsp_snps_pooled.hg38.bed")

BEDTOOLS = "bedtools"
# --------------------------------------------------------------------------------------------

SETS = [("1000GP (MAF ≥ 5%)", TGP_BED),
        ("All shared", SHARED_BED),
        ("Old shared (≥4 Mya)", OLD_SHARED_BED),
        ("TSP", TSP_BED)]


def n_overlap(bed, registry):
    """SNPs in `bed` overlapping any cCRE, and the set size."""
    n = sum(1 for _ in open(bed))
    p = subprocess.run(
        f"cut -f1-3 {bed} | sort -k1,1 -k2,2n | "
        f"{BEDTOOLS} intersect -u -a - -b {registry} | wc -l",
        shell=True, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(p.stderr)
    return int(p.stdout.strip()), n


def cre_classes(bed, registry):
    """Which cCRE classes the overlapping SNPs hit (column 6 of the registry)."""
    p = subprocess.run(
        f"cut -f1-3 {bed} | sort -k1,1 -k2,2n | "
        f"{BEDTOOLS} intersect -wa -wb -a - -b {registry} | cut -f1,2,9 | sort -u | cut -f3",
        shell=True, capture_output=True, text=True)
    return pd.Series(p.stdout.split()).value_counts() if p.stdout.strip() else pd.Series(dtype=int)


def compute(registry, sets=None):
    rows = []
    for name, bed in (sets or SETS):
        obs, n = n_overlap(bed, registry)
        rows.append(dict(snp_set=name, N=n, in_ccre=obs, pct=100 * obs / n))
    return pd.DataFrame(rows)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--registry", default="v4", choices=sorted(REGISTRY),
                    help="SCREEN registry version (default v4, the reported one)")
    ap.add_argument("--tsp-bed", default=TSP_BED, help="alternative TSP set (hg38 BED)")
    ap.add_argument("--out", default=os.path.join(HERE, "snp_set_ccre_overlap.txt"))
    a = ap.parse_args()

    reg = REGISTRY[a.registry]
    print(f"  registry {a.registry}: {sum(1 for _ in open(reg)):,} elements\n"
          f"  TSP set: {os.path.basename(a.tsp_bed)}")

    res = compute(reg, [(n, a.tsp_bed if n == "TSP" else b) for n, b in SETS])
    res.to_csv(a.out, sep="\t", index=False)
    print(res.to_string(index=False))
    print(f"\n  cCRE classes hit by the TSPs:\n"
          f"{cre_classes(a.tsp_bed, reg).to_string()}")
    print(f"\n  wrote {a.out}")
