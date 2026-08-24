#!/usr/bin/env python3
"""Figure 5B: the per-SNP TSP-support track, as a LocusZoom / LocalZoom input.

Panel B is a LocusZoom plot of the MUC7 locus in which the association signal is replaced by
the ARG support for each SNP being a TSP. The plot is made on the LocusZoom site, so this
script only builds the file it is given.

Support is the pooled value: the mean of `tsp_first_or_second_count` over the three seed
chains, i.e. the percentage of all 300 posterior genealogies meeting the topology criteria.
It is written to `neg_log_pvalue`, so the y-axis reads directly as support (%). 
Support is NOT a p-value, and the p-value column is a convenience only.

Three files are written:

    ...allshared.tsv          every shared SNP in the SINGER region (35 over 775 kb)
    ...candidates_window.tsv  candidates only, restricted to the plotted window (primary, plotted in Fig. 5B)
    ...allshared.tsv.gz(.tbi) bgzipped + tabix-indexed, for LocalZoom
    ...candidates_window.tsv.gz(.tbi) bgzipped + tabix-indexed, for LocalZoom
    
LocalZoom column mapping (region queries return no header, so specify by number):
    chromosome 1, position 2, ref 3, alt 4, p-value 5 with "-log10 p-value" ticked.
    Build GRCh37; chromosome names are unprefixed.

Usage:
    python build_panel_b_input.py
"""
import argparse
import os
import subprocess

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "../../.."))

# ---- inputs (edit) -------------------------------------------------------------------------
ANALYSES = "resources"
ARG_RUNS = {s: f"{ANALYSES}/singer_seed{s}" for s in (1, 2, 3)}
# the age >=4 Mya candidate set (stage 3), for the candidates-only file
CANDIDATES = (f"{REPO}/6_TSP_Analyses/snp_set_comparisons/misc_files/"
              "old_shared_snps_860.hg19.chrom_pos.txt")
# --------------------------------------------------------------------------------------------

REGION = "4_70860866_71636175"
WINDOW = (71_150_000, 71_350_000)      # the window of the published panel
COLS = ["chromosome", "position", "ref_allele", "alt_allele", "neg_log_pvalue"]


def pooled_support():
    """Mean tsp_first_or_second_count over the three chains, per shared SNP in the region."""
    per = {}
    for s in sorted(ARG_RUNS):
        f = f"{ARG_RUNS[s]}/singer_trees/{REGION}/{REGION}.all_shared_metrics.grouped.txt"
        t = pd.read_csv(f, sep="\t")
        t = t[t.chrom != "chrom"].copy()          # the file repeats its header
        t["position"] = t.position.astype(int)
        per[f"seed{s}"] = t.set_index("position")["tsp_first_or_second_count"].astype(float)
    m = pd.DataFrame(per)
    m["neg_log_pvalue"] = m.mean(axis=1).round(2)
    return m.reset_index().rename(columns={"index": "position"})


def alleles():
    """REF/ALT for the region, from the merged VCF the ARGs were built from."""
    vcf = f"{ARG_RUNS[1]}/merged_regions/{REGION}.vcf.gz"
    p = subprocess.run(["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\n", vcf],
                       capture_output=True, text=True, check=True)
    rows = [l.split("\t") for l in p.stdout.strip().split("\n")]
    d = pd.DataFrame(rows, columns=["chromosome", "position", "ref_allele", "alt_allele"])
    d["position"] = d.position.astype(int)
    return d.drop_duplicates("position")


def tabix(path):
    """bgzip + tabix, for LocalZoom. 1-based, begin == end, one header line skipped."""
    tmp = path + ".sorted"
    with open(path) as fh, open(tmp, "w") as out:
        head = fh.readline()
        out.write(head)
        out.writelines(sorted(fh, key=lambda l: (l.split("\t")[0], int(l.split("\t")[1]))))
    subprocess.run(f"bgzip -f -c {tmp} > {path}.gz", shell=True, check=True)
    subprocess.run(["tabix", "-f", "-s", "1", "-b", "2", "-e", "2", "-S", "1", f"{path}.gz"],
                   check=True)
    os.remove(tmp)
    return f"{path}.gz"


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out-dir", default=HERE)
    a = ap.parse_args()

    m = pooled_support().merge(alleles(), on="position", how="left")
    m["chromosome"] = m.chromosome.fillna("4")
    m["pvalue"] = np.clip(10.0 ** (-m.neg_log_pvalue), 1e-200, 0.99)

    allshared = os.path.join(a.out_dir, "panelB_support.allshared.tsv")
    m[COLS].to_csv(allshared, sep="\t", index=False)
    print(f"  {len(m)} shared SNPs -> {os.path.basename(allshared)}")

    cand = pd.read_csv(CANDIDATES, sep="\t", header=None, names=["chrom", "pos"], comment="#")
    cand = cand[pd.to_numeric(cand.chrom, errors="coerce").notna()]
    keep = set(cand[cand.chrom.astype(int) == 4].pos.astype(int))
    w = m[m.position.isin(keep) & m.position.between(*WINDOW)]
    wf = os.path.join(a.out_dir, "panelB_support.candidates_window.tsv")
    w[COLS].to_csv(wf, sep="\t", index=False)
    print(f"  {len(w)} candidates in chr4:{WINDOW[0]:,}-{WINDOW[1]:,} -> "
          f"{os.path.basename(wf)}")

    print(f"  tabixed -> {os.path.basename(tabix(allshared))}(.tbi)")
    print(f"\n  SNPs above 80% pooled support:\n"
          f"{m.loc[m.neg_log_pvalue > 80, ['position', 'neg_log_pvalue']].to_string(index=False)}")
