#!/usr/bin/env python3
"""Linkage disequilibrium between the MUC7 TSPs and the archaic MUC7 haplotype.

MUC7 carries a deeply divergent haplotype in African populations that Xu et al. (2017,
Mol Biol Evol 34:2704-2715) attribute to introgression from an unidentified archaic hominin.
The results text states that the TSPs we identify lie close to, but are distinct from and not
in linkage disequilibrium with, that haplotype. This measures it.

The archaic haplotype is the 44 tag SNPs of Xu et al.'s group E (bundled as
xu_2017_data/haplotype_snps.txt), spanning chr4:71,337,964-71,348,157 in hg19. The MUC7 TSPs
sit ~127 kb proximal (the five-TSP region at 71,210,148-71,210,970) and ~38 kb distal (the
single-TSP region at 71,386,175), so they flank the haplotype without overlapping it.

r2 is computed in 1000 Genomes phase 3, within AFR (the populations the haplotype segregates
in) and genome-wide, with PLINK's r2 floor set to 0 so weak pairs are reported rather than
dropped -- the claim is about the absence of LD, so an omitted pair and a low pair must be
distinguishable.

    python muc7_archaic_ld.py

Outputs muc7_archaic_ld.tsv, one row per (TSP, archaic SNP) pair, and prints the maximum r2
per TSP.
"""
import argparse
import os
import subprocess
import tempfile

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))

# ---- inputs (edit) -------------------------------------------------------------------------
# 1000 Genomes phase 3, GRCh37. Not bundled -- see the stage README.
TGP_VCF = ("resources/1000GP/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b."
           "20130502.genotypes.vcf.gz")
TGP_PANEL = "resources/1000GP/integrated_call_samples_v3.20130502.ALL.panel"
ARCHAIC = "xu_2017_data/haplotype_snps.txt"      # Xu et al. 2017 group-E tag SNPs (bundled)
PLINK = "plink"
# --------------------------------------------------------------------------------------------

CHROM = 4
TSPS = [71210148, 71210198, 71210321, 71210336, 71210970, 71386175]
PANELS = {"AFR": ("super_pop", "AFR"), "ALL": (None, None)}


def path(p):
    return p if os.path.isabs(p) else os.path.join(HERE, p)


def samples(panel):
    col, val = PANELS[panel]
    if col is None:
        return None
    p = pd.read_csv(path(TGP_PANEL), sep="\t")
    return sorted(p.loc[p[col] == val, "sample"])


def r2(positions, keep, td):
    """Pairwise r2 among `positions` on CHROM, restricted to `keep` samples."""
    # Restrict to exactly these sites, biallelic SNPs only: extracting the whole span leaves
    # multiallelic records that collide once PLINK derives variant IDs from position.
    vcf = os.path.join(td, "region.vcf.gz")
    reg = os.path.join(td, "sites.txt")
    with open(reg, "w") as fh:
        fh.writelines(f"{CHROM}\t{p}\n" for p in sorted(positions))
    subprocess.run(["bcftools", "view", "-R", reg, "-m2", "-M2", "-v", "snps",
                    "-Oz", "-o", vcf, path(TGP_VCF.format(chrom=CHROM))], check=True,
                   capture_output=True)
    cmd = [PLINK, "--vcf", vcf, "--double-id", "--set-missing-var-ids", "@:#",
           "--r2", "--ld-window", "999999", "--ld-window-kb", "1000",
           "--ld-window-r2", "0",           # report weak pairs; absence must mean absence
           "--out", os.path.join(td, "ld")]
    if keep:
        kf = os.path.join(td, "keep.txt")
        with open(kf, "w") as fh:
            fh.writelines(f"{s}\t{s}\n" for s in keep)
        cmd[6:6] = ["--keep", kf]
    subprocess.run(cmd, check=True, capture_output=True)
    return pd.read_csv(os.path.join(td, "ld.ld"), sep=r"\s+")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default=os.path.join(HERE, "muc7_archaic_ld.tsv"))
    a = ap.parse_args()

    arch = pd.read_csv(path(ARCHAIC), sep="\t")
    apos = sorted(arch["position"].astype(int))
    print(f"  archaic haplotype: {len(apos)} tag SNPs, "
          f"chr{CHROM}:{min(apos):,}-{max(apos):,}")
    print(f"  MUC7 TSPs: {len(TSPS)}, chr{CHROM}:{min(TSPS):,}-{max(TSPS):,}")

    rows = []
    for panel in PANELS:
        keep = samples(panel)
        with tempfile.TemporaryDirectory() as td:
            ld = r2(TSPS + apos, keep, td)
        # keep only TSP x archaic pairs, in either orientation
        t, s = set(TSPS), set(apos)
        m = ld[(ld.BP_A.isin(t) & ld.BP_B.isin(s)) | (ld.BP_B.isin(t) & ld.BP_A.isin(s))].copy()
        m["tsp"] = m.apply(lambda r: r.BP_A if r.BP_A in t else r.BP_B, axis=1)
        m["archaic"] = m.apply(lambda r: r.BP_B if r.BP_A in t else r.BP_A, axis=1)
        m["panel"] = panel
        rows.append(m[["panel", "tsp", "archaic", "R2"]])
        n = f"{len(keep):,}" if keep else "all"
        print(f"\n  {panel} (n = {n} samples): {len(m):,} TSP x archaic pairs")
        for tsp, g in m.groupby("tsp"):
            print(f"    chr{CHROM}:{tsp:<10,}  max r2 = {g.R2.max():.4f}   "
                  f"median {g.R2.median():.4f}")

    out = pd.concat(rows, ignore_index=True)
    out.to_csv(a.out, sep="\t", index=False)
    print(f"\n  wrote {a.out}: {len(out):,} rows")
