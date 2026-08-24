#!/usr/bin/env python3
"""Check that the TSPs within each multi-TSP region are in linkage disequilibrium.

Clustering by distance alone does not establish that a region's TSPs sit on one trans-species
haplotype rather than being independent nearby polymorphisms. This measures r^2 between every
pair of TSPs in a region, in the 1000 Genomes phase 3 panel, within three sample sets:

    AFR  all African samples          (n = 661)
    EUR  all European samples         (n = 503)
    YRI  Yoruba in Ibadan only        (n = 108)

AFR and EUR are the panels reported; YRI is included because it is the population whose Relate
ages define the TSP set, and it is the smallest, so it is the weakest test.

For each region and panel: the minimum and mean pairwise r^2 across all TSP pairs, and whether
the minimum reaches 0.8 (one tight block). A region where the minimum is low but pairs are
individually high contains more than one block; IGFBP7 is the example in the published set.

Only multi-TSP regions are tested. Singletons have no pairs.

Coordinates are hg19/b37 throughout, matching the phase 3 panel.

Inputs:  tsp_clusters_10kb.pooled.tsv and tsps.pooled.hg19.txt from aggregate_seeds.py;
         the 1000 Genomes phase 3 per-chromosome VCFs and the sample panel file;
         optionally ../3_GeneAssignment/region_gene_assignments.tsv, for the gene label only.
Output:  within_region_LD.txt

The inputs below are the published defaults. They can be overridden on the command line
(--clusters/--tsps/--out/--genes/--tgp-vcf/--tgp-panel) to run the same test on an
alternative TSP set, which is how the Ne-sensitivity comparison uses it.
"""
import argparse
import os
import subprocess
import tempfile

import numpy as np
import pandas as pd

# ---- inputs (edit) -------------------------------------------------------------------------
CLUSTERS = "tsp_clusters_10kb.pooled.tsv"
TSPS = "tsps.pooled.hg19.txt"

# 1000 Genomes phase 3, GRCh37. {chrom} is substituted.
TGP_VCF = ("resources/1000GP/"
           "ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz")
# The phase 3 sample panel (sample / pop / super_pop / gender), from the same release:
# https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
TGP_PANEL = "resources/1000GP/integrated_call_samples_v3.20130502.ALL.panel"

# Gene labels, for readability only; regenerate after ../3_GeneAssignment/ has run.
GENE_ASSIGNMENTS = "../3_GeneAssignment/region_gene_assignments.pooled.tsv"

OUT = "within_region_LD.txt"
R2_THRESHOLD = 0.8
BCFTOOLS = "bcftools"
PLINK = "plink"
# --------------------------------------------------------------------------------------------

# panel -> the three sample sets. AFR/EUR are super-populations, YRI a single population.
PANELS = {"AFR": ("super_pop", "AFR"), "EUR": ("super_pop", "EUR"), "YRI": ("pop", "YRI")}

HERE = os.path.dirname(os.path.abspath(__file__))


def path(p):
    return p if os.path.isabs(p) else os.path.join(HERE, p)


def sample_sets():
    panel = pd.read_csv(path(TGP_PANEL), sep="\t")
    return {name: sorted(panel.loc[panel[col] == value, "sample"])
            for name, (col, value) in PANELS.items()}


def pairwise_r2(chrom, positions, keep_samples, tmpdir):
    """All pairwise r^2 among `positions` on `chrom`, within `keep_samples`.

    Returns (r2 values, number of the requested positions actually present in the panel).
    """
    prefix = os.path.join(tmpdir, "region")
    sites = f"{prefix}.sites"
    with open(sites, "w") as f:
        for p in positions:
            f.write(f"{chrom}\t{p}\n")
    keep = f"{prefix}.keep"
    with open(keep, "w") as f:
        for s in keep_samples:
            f.write(f"{s}\t{s}\n")     # --double-id sets FID = IID = sample

    vcf = f"{prefix}.vcf.gz"
    subprocess.run([BCFTOOLS, "view", "-m2", "-M2", "-v", "snps", "-R", sites,
                    "-Oz", "-o", vcf, path(TGP_VCF.format(chrom=chrom))], check=True)

    # --ld-window-r2 0 so weakly correlated pairs are reported too; without it PLINK's default
    # 0.2 filter would silently drop exactly the pairs this check exists to find.
    run = subprocess.run(
        [PLINK, "--vcf", vcf, "--double-id", "--keep", keep,
         "--set-missing-var-ids", "@:#", "--r2",
         "--ld-window", "99999", "--ld-window-kb", "1000", "--ld-window-r2", "0",
         "--out", prefix],
        capture_output=True, text=True)
    if run.returncode != 0:
        raise RuntimeError(run.stdout + run.stderr)

    n_in_vcf = int(subprocess.run([BCFTOOLS, "view", "-H", vcf], capture_output=True,
                                  text=True, check=True).stdout.count("\n"))
    ld_file = f"{prefix}.ld"
    if not os.path.exists(ld_file):
        return np.array([]), n_in_vcf
    ld = pd.read_csv(ld_file, sep=r"\s+")
    return ld["R2"].to_numpy(dtype=float), n_in_vcf


def parse_args():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--clusters", default=CLUSTERS, help="per-region table from cluster_tsps.py")
    p.add_argument("--tsps", default=TSPS, help="per-TSP table with the `region` column")
    p.add_argument("--out", default=OUT)
    p.add_argument("--genes", default=GENE_ASSIGNMENTS,
                   help="region_gene_assignments.tsv; used for the label column only")
    p.add_argument("--tgp-vcf", default=TGP_VCF, help="1000GP phase 3 VCF, with {chrom}")
    p.add_argument("--tgp-panel", default=TGP_PANEL)
    p.add_argument("--quiet", action="store_true", help="suppress the per-row table on stdout")
    p.add_argument("--tsv", help="also write a tab-separated copy here. The default output is "
                                 "an aligned fixed-width table, which is not safely parseable "
                                 "when a gene label is blank.")
    return p.parse_args()


def main():
    global CLUSTERS, TSPS, OUT, GENE_ASSIGNMENTS, TGP_VCF, TGP_PANEL
    a = parse_args()
    CLUSTERS, TSPS, OUT = a.clusters, a.tsps, a.out
    GENE_ASSIGNMENTS, TGP_VCF, TGP_PANEL = a.genes, a.tgp_vcf, a.tgp_panel

    clusters = pd.read_csv(path(CLUSTERS), sep="\t")
    tsps = pd.read_csv(path(TSPS), sep="\t", low_memory=False)
    samples = sample_sets()
    for name, s in samples.items():
        print(f"  {name}: {len(s)} samples")

    genes = {}
    if os.path.exists(path(GENE_ASSIGNMENTS)):
        ga = pd.read_csv(path(GENE_ASSIGNMENTS), sep="\t")
        genes = dict(zip(ga["region"], ga["assigned_gene"]))
    else:
        print(f"  note: {GENE_ASSIGNMENTS} not found, gene column left blank")

    multi = clusters[clusters["n_tsps"] > 1].sort_values("n_tsps", ascending=False)
    print(f"testing {len(multi)} multi-TSP regions "
          f"({len(clusters) - len(multi)} singletons have no pairs)\n")

    rows = []
    for c in multi.itertuples():
        region = f"{c.chrom}_{c.start}_{c.end}"
        positions = sorted(tsps.loc[tsps["region"] == region, "POS"])
        for panel_name, keep in samples.items():
            with tempfile.TemporaryDirectory() as tmpdir:
                r2, n_in_vcf = pairwise_r2(c.chrom, positions, keep, tmpdir)
            rows.append({
                "region": region,
                "gene": genes.get(region, ""),
                "nSNP": len(positions),
                "nInVCF": n_in_vcf,
                "pop": panel_name,
                "minR2": round(float(r2.min()), 3) if len(r2) else np.nan,
                "meanR2": round(float(r2.mean()), 3) if len(r2) else np.nan,
                ">=0.8?": ("yes" if len(r2) and r2.min() >= R2_THRESHOLD else "NO"),
            })

    out = pd.DataFrame(rows)
    with open(path(OUT), "w") as f:
        f.write(out.to_string(index=False) + "\n")
    if a.tsv:
        out.to_csv(path(a.tsv), sep="\t", index=False)
    if not a.quiet:
        print(out.to_string(index=False))
    print(f"\nwrote {OUT}" + (f" and {a.tsv}" if a.tsv else ""))

    failed = out[out[">=0.8?"] == "NO"]
    if len(failed):
        print(f"\n{len(failed)} region x panel combinations below r2 = {R2_THRESHOLD}:")
        print(failed[["region", "gene", "pop", "minR2", "meanR2"]].to_string(index=False))


if __name__ == "__main__":
    main()
