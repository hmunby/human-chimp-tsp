#!/usr/bin/env python3
"""phyloP conservation: descriptive statistics and tests across the SNP sets.

Two-sided Mann-Whitney (with a common-language effect size and Cliff's delta) plus a
two-sample Kolmogorov-Smirnov test for each pair of interest, Benjamini-Hochberg corrected
across the Mann-Whitney tests. Reported with and without the MHC so the filter can be
checked against the rest of the paper.

Sets (hg19): chimpanzee AC>=3, all shared, old shared (>=4 Mya), TSP, 1000GP MAF>=5%.
Tracks:      phylop/       36-eutherian-mammal EPO, human EXCLUDED (main figure)
             phylop_46way/ UCSC 46-way placental, human INCLUDED (comparison)

Inputs:  phylop*/ per-chromosome score tables from ./Snakefile; the SNP-set position lists
         in misc_files/.
Output:  the table printed to stdout (redirect to phylop_stats_output.txt).
"""
import os

import numpy as np
import pandas as pd
from scipy.stats import ks_2samp, mannwhitneyu

HERE = os.path.dirname(os.path.abspath(__file__))

# ---- inputs (edit) -------------------------------------------------------------------------
OLD_SHARED_FILE = os.path.join(HERE, "misc_files/old_shared_snps_860.hg19.chrom_pos.txt")
TSP_SNP_FILE = os.path.join(HERE, "misc_files/tsp_snps_pooled.hg19.txt")
MHC_CHROM, MHC_START, MHC_END = 6, 28477797, 33448354   # MHC, hg19 (these tables are hg19-keyed)
# --------------------------------------------------------------------------------------------

LABEL = {"chimp": "Chimpanzee AC>=3", "allshared": "All shared", "old4": ">4 Mya shared",
         "tsp": "TSP", "tgp": "1000GP MAF>=5%"}

PAIRS = [("tsp", "tgp"), ("tsp", "allshared"), ("tsp", "old4"), ("tsp", "chimp"),
         ("allshared", "tgp"), ("chimp", "tgp"), ("old4", "tgp")]


def _load_chrom_pos(path):
    df = pd.read_csv(path, sep="\t", low_memory=False)
    df = df.rename(columns={c: c.upper() for c in df.columns})
    df["CHROM"] = pd.to_numeric(df["CHROM"], errors="coerce")
    df["POS"] = pd.to_numeric(df["POS"], errors="coerce")
    return df.dropna(subset=["CHROM", "POS"]).astype({"CHROM": int, "POS": int})


def load_tracks(phylop_dir, drop_mhc, tsp_file=None, old_shared_file=None):
    """The five phyloP tracks.

    tsp_file/old_shared_file default to the published sets; passing others scores an
    alternative TSP set through identical code (the Ne-sensitivity comparison does that).
    """
    tsp_file = tsp_file or TSP_SNP_FILE
    old_shared_file = old_shared_file or OLD_SHARED_FILE
    d = os.path.join(HERE, phylop_dir)
    shared = {c: pd.read_csv(f"{d}/shared_snps_phylop_chr{c}.txt", sep="\t") for c in range(1, 23)}

    def mask_mhc(df):
        if not drop_mhc:
            return df
        bad = (df["CHROM"] == MHC_CHROM) & df["POS"].between(MHC_START, MHC_END)
        return df[~bad]

    def flat(prefix):
        parts = [mask_mhc(pd.read_csv(f"{d}/{prefix}_chr{c}.txt", sep="\t")) for c in range(1, 23)]
        return pd.concat(parts)["PHYLOP_SCORE"].to_numpy(dtype=float)

    def subset(chrom_pos):
        out = []
        for c in range(1, 23):
            sub = chrom_pos[chrom_pos["CHROM"] == c]
            if len(sub):
                m = pd.merge(sub[["CHROM", "POS"]], shared[c], on=["CHROM", "POS"], how="inner")
                out.append(mask_mhc(m))
        return pd.concat(out)["PHYLOP_SCORE"].to_numpy(dtype=float) if out else np.array([])

    t = {"chimp": flat("chimp_snps_phylop"),
         "allshared": pd.concat([mask_mhc(shared[c]) for c in range(1, 23)])["PHYLOP_SCORE"].to_numpy(float),
         "old4": subset(_load_chrom_pos(old_shared_file)),
         "tsp": subset(_load_chrom_pos(tsp_file)),
         "tgp": flat("tgp_maf5_snps_phylop")}
    return {k: v[~np.isnan(v)] for k, v in t.items()}


def describe(tracks):
    print(f"  {'set':18} {'n':>9} {'mean':>8} {'median':>8} {'IQR':>18} {'% > 0':>7}")
    for k, v in tracks.items():
        q1, q3 = np.percentile(v, [25, 75])
        print(f"  {LABEL[k]:18} {len(v):>9,} {np.mean(v):>+8.3f} {np.median(v):>+8.3f}"
              f"  [{q1:+.2f}, {q3:+.2f}]{'':>4} {100 * (v > 0).mean():>6.1f}")


def test(tracks):
    rows = []
    for a, b in PAIRS:
        x, y = tracks[a], tracks[b]
        U, p_mw = mannwhitneyu(x, y, alternative="two-sided")
        auc = U / (len(x) * len(y))     # P(x > y), common-language effect size
        D, p_ks = ks_2samp(x, y)
        rows.append(dict(a=LABEL[a], b=LABEL[b], n_a=len(x), n_b=len(y),
                         d_median=np.median(x) - np.median(y), auc=auc, cliff=2 * auc - 1,
                         p_mw=p_mw, D=D, p_ks=p_ks))
    df = pd.DataFrame(rows)

    # Benjamini-Hochberg across the Mann-Whitney tests
    p = df["p_mw"].to_numpy()
    order = np.argsort(p)
    m = len(df)
    q = np.empty(m)
    running = 1.0
    for rank, i in enumerate(order[::-1]):
        running = min(running, p[i] * m / (m - rank))
        q[i] = running
    df["q_mw"] = q

    print(f"\n  {'comparison':40} {'Δmedian':>8} {'P(a>b)':>7} {'Cliff δ':>8} {'MWU p':>10}"
          f" {'MWU q':>10} {'KS D':>6} {'KS p':>10}")
    for _, r in df.iterrows():
        print(f"  {r['a'] + ' vs ' + r['b']:40} {r['d_median']:>+8.3f} {r['auc']:>7.3f}"
              f" {r['cliff']:>+8.3f} {r['p_mw']:>10.2e} {r['q_mw']:>10.2e} {r['D']:>6.3f}"
              f" {r['p_ks']:>10.2e}")
    return df


if __name__ == "__main__":
    import argparse
    _ap = argparse.ArgumentParser(description=__doc__.split("\n")[0] if __doc__ else None)
    _ap.add_argument("--tsp-file", help="alternative TSP set; defaults to the reported pooled set")
    _ap.add_argument("--old-shared-file", help="alternative old-shared set")
    _a = _ap.parse_args()

    for phylop_dir, version in [("phylop", "36-mammal EPO, HUMAN EXCLUDED (main)"),
                                ("phylop_46way", "UCSC 46-way, human included")]:
        for drop_mhc in (False, True):
            print("\n" + "=" * 100)
            print(f"{version}   |   MHC {'EXCLUDED' if drop_mhc else 'included'}")
            print("=" * 100)
            tr = load_tracks(phylop_dir, drop_mhc, _a.tsp_file, _a.old_shared_file)
            describe(tr)
            test(tr)
