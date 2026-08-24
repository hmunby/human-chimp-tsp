#!/usr/bin/env python3
"""Combine the three Ne=50,000 MCMC chains into the reported TSP set.

Each chain (../../3_SINGER, runs ne50000_seed1/2/3) scores every shared SNP with a
support count out of 100 posterior trees. This applies the TSP test at the original
candidate SNPs and clusters the survivors into the reported regions:

    candidate SNP  ->  support rule  ->  single-linkage cluster at 10 kb  ->  regions

Nothing else is imposed. In particular the age criterion is NOT re-applied: it is what
selected the candidates in the first place (shared SNPs with a YRI midpoint age over
4 Mya), so filtering on age again after the ARG test applies the same cut twice.

## The rules

  pooled     REPORTED. Mean support across the three chains, thresholded once. Because
             each chain contributes exactly 100 trees, this is identical to pooling all
             300 posterior trees and asking whether more than 80% of them show the
             trans-species configuration -- i.e. treating the three chains as the one
             longer chain they are.
  strict     Sensitivity check. Passes independently in all three chains.
  majority   Passes in at least two of three.
  any        Passes in at least one.

`pooled` is the reported rule; `strict` is written as the sensitivity set. Both, and the
other two, are produced by one invocation, so any of them can be passed downstream.

Note the threshold convention: the count-based rules are sensitive to `>` versus `>=`,
because a single chain can land on exactly 80 whereas a mean of three rarely does.
`--inclusive` switches it.

Inputs:  ../../3_SINGER/singer/<run>/*/*.all_shared_metrics.grouped.txt
         the candidate SNP list from ../../3_SINGER/misc_files/
Outputs: tsps.<rule>.hg19.txt        the TSP table, `region` = reported cluster id
         tsps.<rule>.hg38.bed        hg38 positions (annotation; nothing is filtered on it)
         tsp_clusters_10kb.<rule>.tsv    one row per reported region
         seed_support.tsv            per candidate: support in each chain, and every verdict
"""
import argparse
import glob
import os
import subprocess
import tempfile

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))

# ---- inputs (edit) -------------------------------------------------------------------------
# One entry per chain: a run name from ../../3_SINGER/Snakefile's RUNS dict. All three are
# produced by the same workflow, differing only in the MCMC seed.
SINGER = "../../3_SINGER/singer"
SEED_RUNS = {
    "seed1": f"{SINGER}/ne50000_seed1",
    "seed2": f"{SINGER}/ne50000_seed2",
    "seed3": f"{SINGER}/ne50000_seed3",
}
CANDIDATES = "../../3_SINGER/misc_files/merged_candidate_snps.chrom_pos.txt"
REGION_IDS = "../../3_SINGER/misc_files/non_overlapping_region_250kb_ids.txt"
EXCLUDED = "../../3_SINGER/misc_files/excluded_regions.txt"
# Shared-SNP annotation table (alleles, AF, BMAP, CpG, mutation rate, Relate ages). Large;
# not bundled -- see the stage README.
SNP_INFO = "resources/all_shared_snps_ann.txt"
CHAIN = "resources/b37tohg38.nochr.over.chain"

THRESHOLD = 80         # support strictly greater than this, of 100 posterior trees
INCLUSIVE = False      # True applies >= instead of >
CLUSTER_GAP = 10_000   # single-linkage gap, bp
GEN_TIME = 29
# --------------------------------------------------------------------------------------------

GENE_ASSIGNMENT_COLS = ["CHROM_HG38", "POSITION_HG38", "REF", "ALT", "BMAP", "Annotation",
                        "Gene_Name", "region", "CHROM", "POS", "MID_AGE_YRI_YEARS",
                        "tsp_first_or_second_count"]


def path(p):
    return p if os.path.isabs(p) else os.path.join(HERE, p)


def valid_regions():
    target = {l.strip() for l in open(path(REGION_IDS)) if l.strip()}
    excl = set()
    if os.path.exists(path(EXCLUDED)):
        excl = {l.split()[0] for l in open(path(EXCLUDED)) if l.strip() and not l.startswith("#")}
    return target - excl


def completed(run_dir, regions):
    """Regions this chain has finished (grouped metrics present)."""
    return {r for r in regions
            if os.path.exists(f"{path(run_dir)}/{r}/{r}.all_shared_metrics.grouped.txt")}


def candidates():
    """The original candidate SNPs (autosomes). This IS the age criterion."""
    out = set()
    for line in open(path(CANDIDATES)):
        p = line.split()
        if len(p) >= 2 and p[0] != "X":
            out.add((int(p[0]), int(p[1])))
    return out


def support_table():
    regions = valid_regions()
    comp = {k: completed(v, regions) for k, v in SEED_RUNS.items()}
    common = set.intersection(*comp.values())
    print("regions completed: " + "  ".join(f"{k}={len(v)}" for k, v in comp.items()))
    print(f"common to all three chains: {len(common)}")

    cand = candidates()
    frames = {}
    for k, d in SEED_RUNS.items():
        rows = []
        for r in sorted(common):
            f = f"{path(d)}/{r}/{r}.all_shared_metrics.grouped.txt"
            try:
                m = pd.read_csv(f, sep="\t")
            except Exception:
                continue
            if len(m):
                rows.append(m[["chrom", "position", "tsp_first_or_second_count"]])
        m = pd.concat(rows, ignore_index=True).drop_duplicates(subset=["chrom", "position"])
        m = m[[(int(c), int(p)) in cand for c, p in zip(m["chrom"], m["position"])]]
        frames[k] = m.rename(columns={"tsp_first_or_second_count": k})

    d = frames["seed1"]
    for k in ("seed2", "seed3"):
        d = d.merge(frames[k], on=["chrom", "position"], how="inner")
    d["mean_support"] = d[["seed1", "seed2", "seed3"]].mean(axis=1)
    print(f"candidates tested in all three chains: {len(d)}")
    return d


def verdicts(d):
    hit = (lambda x: x >= THRESHOLD) if INCLUSIVE else (lambda x: x > THRESHOLD)
    p = {k: hit(d[k]) for k in ("seed1", "seed2", "seed3")}
    n = sum(p[k].astype(int) for k in p)
    return {**p, "any": n >= 1, "majority": n >= 2, "strict": n == 3,
            "pooled": hit(d["mean_support"])}


def cluster(t):
    """Single-linkage at CLUSTER_GAP; region id is chrom_start_end over its own TSPs."""
    t = t.sort_values(["CHROM", "POS"]).reset_index(drop=True)
    ids, nxt = [], 0
    for _, g in t.groupby("CHROM", sort=True):
        prev = None
        for pos in g["POS"]:
            if prev is not None and pos - prev > CLUSTER_GAP:
                nxt += 1
            ids.append(nxt)
            prev = pos
        nxt += 1
    t["_c"] = ids
    cl = t.groupby("_c").agg(chrom=("CHROM", "first"), start=("POS", "min"),
                             end=("POS", "max"), n_tsps=("POS", "size"),
                             min_support=("tsp_first_or_second_count", "min"),
                             max_support=("tsp_first_or_second_count", "max")).reset_index(drop=True)
    cl["region"] = (cl["chrom"].astype(str) + "_" + cl["start"].astype(str)
                    + "_" + cl["end"].astype(str))
    cl["span_bp"] = cl["end"] - cl["start"]
    t["region"] = t["_c"].map(t.groupby("_c").apply(
        lambda g: f"{g['CHROM'].iloc[0]}_{g['POS'].min()}_{g['POS'].max()}",
        include_groups=False))
    cols = ["region", "chrom", "start", "end", "span_bp", "n_tsps", "min_support", "max_support"]
    return t.drop(columns="_c"), cl[cols]


def liftover(t):
    """Annotate hg38 coordinates. Nothing is filtered on this."""
    with tempfile.TemporaryDirectory() as td:
        h19, h38, un = (os.path.join(td, x) for x in ("a.bed", "b.bed", "u.bed"))
        b = t[["CHROM", "POS"]].copy()
        b["start"] = b["POS"] - 1
        b["name"] = b["CHROM"].astype(str) + "_" + b["POS"].astype(str)
        b[["CHROM", "start", "POS", "name"]].to_csv(h19, sep="\t", header=False, index=False)
        subprocess.run(["liftOver", h19, path(CHAIN), h38, un], check=False,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if os.path.exists(h38) and os.path.getsize(h38):
            lift = pd.read_csv(h38, sep="\t", header=None, names=["c", "s", "e", "name"])
            lift["CHROM_HG38"] = lift["c"].astype(str).str.replace("chr", "", regex=False)
            lift["POSITION_HG38"] = lift["e"]
            lmap = lift.set_index("name")[["CHROM_HG38", "POSITION_HG38"]]
        else:
            lmap = pd.DataFrame(columns=["CHROM_HG38", "POSITION_HG38"])
    t = t.copy()
    t["name"] = t["CHROM"].astype(str) + "_" + t["POS"].astype(str)
    return t.join(lmap, on="name").drop(columns="name")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--rules", nargs="*", default=["pooled", "strict", "majority", "any",
                                                   "seed1", "seed2", "seed3"])
    ap.add_argument("--threshold", type=int, default=THRESHOLD)
    ap.add_argument("--inclusive", action="store_true", help="apply >= instead of >")
    ap.add_argument("--out-dir", default=HERE)
    a = ap.parse_args()
    globals()["THRESHOLD"], globals()["INCLUSIVE"] = a.threshold, a.inclusive
    os.makedirs(a.out_dir, exist_ok=True)

    d = support_table()
    sel = verdicts(d)
    for k, v in sel.items():
        d[f"is_{k}"] = v
    d.to_csv(os.path.join(a.out_dir, "seed_support.tsv"), sep="\t", index=False)

    info = pd.read_csv(path(SNP_INFO), sep="\t", low_memory=False)
    for rule in a.rules:
        sub = d[sel[rule]].copy()
        sub["tsp_first_or_second_count"] = (
            sub[rule] if rule.startswith("seed") else sub["mean_support"]).round().astype(int)
        t = info.merge(sub, left_on=["CHROM", "POS"], right_on=["chrom", "position"], how="right")
        lo = pd.to_numeric(t["LOWER_AGE_YRI"], errors="coerce")
        up = pd.to_numeric(t["UPPER_AGE_YRI"], errors="coerce")
        t["MID_AGE_YRI_YEARS"] = (lo + up) / 2 * GEN_TIME     # descriptive only
        t = liftover(t)
        t, cl = cluster(t)

        t.to_csv(f"{a.out_dir}/tsps.{rule}.hg19.txt", sep="\t", index=False)
        cl.to_csv(f"{a.out_dir}/tsp_clusters_10kb.{rule}.tsv", sep="\t", index=False)
        have = [c for c in GENE_ASSIGNMENT_COLS if c in t.columns]
        t[have].to_csv(f"{a.out_dir}/tsps.{rule}.for_gene_assignment.txt", sep="\t", index=False)
        b = t.dropna(subset=["POSITION_HG38"]).copy()
        if len(b):
            b["start"] = b["POSITION_HG38"].astype(int) - 1
            b["chrom"] = "chr" + b["CHROM_HG38"].astype(str)
            b["_c"] = pd.to_numeric(b["CHROM_HG38"], errors="coerce")
            b.sort_values(["_c", "start"])[["chrom", "start", "POSITION_HG38"]].to_csv(
                f"{a.out_dir}/tsps.{rule}.hg38.bed", sep="\t", header=False, index=False)
        print(f"  {rule:<9} {len(t):>4} TSPs  {len(cl):>3} regions "
              f"({int((cl['n_tsps'] > 1).sum())} multi-TSP)")


if __name__ == "__main__":
    main()
