#!/usr/bin/env python3
"""Call TSPs from a SINGLE SINGER chain and cluster them into reported regions.

The reported set comes from three chains combined (`aggregate_seeds.py`). This is the
one-chain equivalent, used for the sensitivity runs at Ne = 20,000 and Ne = 100,000, and
it applies exactly the same rule so those sets are comparable to the reported one:

    candidate SNP  ->  support > 80 of 100 posterior trees  ->  cluster at 10 kb

Nothing else is imposed. In particular the age criterion is NOT re-applied: it is what
selected the candidates in the first place (shared SNPs with a YRI midpoint age over
4 Mya), so filtering on age again after the ARG test applies the same cut twice.

Inputs:  ../../3_SINGER/results/<run>/ -- per-region grouped metrics for one run
         the candidate SNP list from ../../1_RecallRegions/
Outputs: tsps.<tag>.hg19.txt, tsp_clusters_10kb.<tag>.tsv, tsps.<tag>.for_gene_assignment.txt

Usage:
    python cluster_tsps.py --run ../../3_SINGER --tag ne20k
    python cluster_tsps.py --run ../../3_SINGER_Ne100k --tag ne100k --threshold 80
"""
import argparse
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
# stage A/B helpers are shared with the multi-chain path so the two cannot drift
from aggregate_seeds import (CLUSTER_GAP, GEN_TIME, GENE_ASSIGNMENT_COLS, SNP_INFO,  # noqa: E402
                             candidates, cluster, liftover, path, valid_regions)


def load_run(run_dir, regions):
    """Per-SNP support from one chain's grouped metrics."""
    rows = []
    for r in sorted(regions):
        f = f"{path(run_dir)}/singer_trees/{r}/{r}.all_shared_metrics.grouped.txt"
        if not os.path.exists(f):
            continue
        try:
            m = pd.read_csv(f, sep="\t")
        except Exception:
            continue
        if len(m):
            rows.append(m[["chrom", "position", "tsp_first_or_second_count"]])
    if not rows:
        raise SystemExit(f"no grouped metrics found under {run_dir}")
    return pd.concat(rows, ignore_index=True).drop_duplicates(subset=["chrom", "position"])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run", required=True, help="a completed SINGER run directory")
    ap.add_argument("--tag", required=True, help="name for the outputs, e.g. ne20k")
    ap.add_argument("--threshold", type=int, default=80)
    ap.add_argument("--inclusive", action="store_true", help="apply >= instead of >")
    ap.add_argument("--out-dir", default=HERE)
    a = ap.parse_args()

    regions = valid_regions()
    d = load_run(a.run, regions)
    cand = candidates()
    d = d[[(int(c), int(p)) in cand for c, p in zip(d["chrom"], d["position"])]]
    print(f"candidates scored in {a.run}: {len(d)}")

    sup = d["tsp_first_or_second_count"]
    d = d[sup >= a.threshold if a.inclusive else sup > a.threshold].copy()
    print(f"passing support {'>=' if a.inclusive else '>'}{a.threshold}: {len(d)}")

    info = pd.read_csv(path(SNP_INFO), sep="\t", low_memory=False)
    t = info.merge(d, left_on=["CHROM", "POS"], right_on=["chrom", "position"], how="right")
    lo = pd.to_numeric(t["LOWER_AGE_YRI"], errors="coerce")
    up = pd.to_numeric(t["UPPER_AGE_YRI"], errors="coerce")
    t["MID_AGE_YRI_YEARS"] = (lo + up) / 2 * GEN_TIME     # descriptive only
    t = liftover(t)
    t, cl = cluster(t)

    os.makedirs(a.out_dir, exist_ok=True)
    t.to_csv(f"{a.out_dir}/tsps.{a.tag}.hg19.txt", sep="\t", index=False)
    cl.to_csv(f"{a.out_dir}/tsp_clusters_10kb.{a.tag}.tsv", sep="\t", index=False)
    have = [c for c in GENE_ASSIGNMENT_COLS if c in t.columns]
    t[have].to_csv(f"{a.out_dir}/tsps.{a.tag}.for_gene_assignment.txt", sep="\t", index=False)
    print(f"  {a.tag}: {len(t)} TSPs in {len(cl)} regions "
          f"({int((cl['n_tsps'] > 1).sum())} multi-TSP)")


if __name__ == "__main__":
    main()
