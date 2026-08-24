#!/usr/bin/env python3
"""
SINGER analogue of main Figure 1A: percentage of SNPs older than age, for
shared SNPs vs all common 1000GP SNPs (MAF>=5%), using SINGER posterior ages
(integrating over the 100 posterior ARG samples per SNP).

Source of 100-sample posterior ages (hg38-keyed, generations):
  resources/african_chr{N}_ages.pickle
    dict: pos(hg38,int) -> list of 100 (lower_gen, upper_gen) tuples

SNP sets (hg38 coordinates):
  shared      : the human-chimp shared SNPs (hg38); those without a SINGER age drop out at the
                pickle lookup, so the denominator is the number actually dated
  all-common  : all_snps.hg38.singer.txt rows with SINGER_AFR==1 (1000GP MAF>=5% with a SINGER age)

Model: under neutrality a mutation is uniform along its branch, so for one posterior
sample tuple (a0,a1) the prob the allele is older than t is clip((a1-t)/(a1-a0),0,1),
i.e. the survival function of Uniform(a0,a1). Averaging over the 100 samples gives each
SNP mass 1. Summed over SNPs and divided by N gives the expected FRACTION older than t.

We accumulate this survival curve on a generation grid with a difference-array over the
uniform densities (each tuple weight 1/100), so the whole thing is O(#tuples + #grid).
"""
import pickle, numpy as np, os, sys

AGE_DIR   = "resources/ages_all"
SHARED    = "resources/shared_snps_chrom_pos_hg38.sorted.bed"
ALLCOMMON = "resources/all_snps.hg38.singer.txt"
OUT       = "outputs"

DT       = 1000                 # grid spacing, generations
MAXGEN   = 2_000_000            # covers the full SINGER age range (~1.96M gen)
NSAMP    = 100
GRID     = np.arange(0, MAXGEN + DT, DT, dtype=np.float64)
NBIN     = len(GRID)

def load_shared_positions():
    """hg38 positions of the shared SNPs, per chrom -> set(int).

    A BED of CHROM/START/END; END is the 1-based position. Positions with no SINGER age are
    counted as missing at the pickle lookup and do not enter the denominator.
    """
    d = {c: set() for c in range(1, 23)}
    with open(SHARED) as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 3:
                continue
            ch = c[0][3:] if c[0].startswith("chr") else c[0]
            try:
                chrom = int(ch)
            except ValueError:
                continue
            if 1 <= chrom <= 22:
                d[chrom].add(int(c[2]))
    return d

def load_common_positions(maf_min=0.05):
    """hg38 positions of 1000GP MAF>=maf_min SNPs with a SINGER age, per chrom -> set(int)."""
    d = {c: set() for c in range(1, 23)}
    with open(ALLCOMMON) as f:
        h = f.readline().rstrip("\n").split("\t"); ci = {n: i for i, n in enumerate(h)}
        af_i, sg_i, ch_i, ps_i = ci["AF"], ci["SINGER_AFR"], ci["CHROM"], ci["POS"]
        for line in f:
            c = line.rstrip("\n").split("\t")
            if c[sg_i] != "1":
                continue
            ch = c[ch_i]
            if not ch.startswith("chr"): continue
            try: chrom = int(ch[3:])
            except ValueError: continue
            if not (1 <= chrom <= 22): continue
            try: af = float(c[af_i])
            except ValueError: continue
            if min(af, 1 - af) < maf_min:
                continue
            d[chrom].add(int(c[ps_i]))
    return d

def _prefix_at(sorted_key, prefix, grid):
    """For each t in grid, return prefix-sum over items whose key < t. Exact."""
    k = np.searchsorted(sorted_key, grid, side='left')      # count of key < t
    out = np.zeros(len(grid))
    nz = k > 0
    out[nz] = prefix[k[nz] - 1]
    return out

def accumulate_chrom_exact(survival, a0, a1):
    """Add one chromosome's posterior tuples to the (unnormalized) survival curve
    S(t) = sum_tuples clip((a1-t)/(a1-a0), 0, 1) / NSAMP, evaluated exactly on GRID.

    Decompose per tuple: ramp=1 if a0>=t (fully older); ramp=(a1-t)/(a1-a0) if
    a0<t<=a1; ramp=0 if a1<t. Since a1>a0, {a1<t} is a subset of {a0<t}, so
       sum_{a0<t<=a1} f = sum_{a0<t} f - sum_{a1<t} f
    with f = a1*w - t*w and w = 1/(a1-a0). All sums are exact cumulative sums."""
    # "fully older" count uses ALL tuples' a0 (good + zero-length treated as point mass at a0)
    alive = len(a0) - np.searchsorted(np.sort(a0), GRID, side='left')

    good = a1 > a0
    active = np.zeros(NBIN)
    if good.any():
        g0 = a0[good]; g1 = a1[good]
        w = 1.0 / (g1 - g0); aw = g1 * w
        o0 = np.argsort(g0)
        sumA0_aw = _prefix_at(g0[o0], np.cumsum(aw[o0]), GRID)
        sumA0_w  = _prefix_at(g0[o0], np.cumsum(w[o0]),  GRID)
        o1 = np.argsort(g1)
        sumA1_aw = _prefix_at(g1[o1], np.cumsum(aw[o1]), GRID)
        sumA1_w  = _prefix_at(g1[o1], np.cumsum(w[o1]),  GRID)
        active = (sumA0_aw - sumA1_aw) - GRID * (sumA0_w - sumA1_w)

    survival += (alive + active) / NSAMP

def run_set(name, positions_by_chrom, maf_note):
    survival = np.zeros(NBIN, dtype=np.float64)
    n_snps = 0; n_missing = 0
    for chrom in range(1, 23):
        pos_set = positions_by_chrom[chrom]
        if not pos_set: continue
        ages = pickle.load(open(f"{AGE_DIR}/african_chr{chrom}_ages.pickle", "rb"))
        a0_list = []; a1_list = []; got = 0; miss = 0
        for p in pos_set:
            a = ages.get(p)
            if a is None:
                miss += 1; continue
            arr = np.asarray(a, dtype=np.float64)   # (100,2)
            a0_list.append(arr[:, 0]); a1_list.append(arr[:, 1])
            got += 1
        if a0_list:
            accumulate_chrom_exact(survival, np.concatenate(a0_list), np.concatenate(a1_list))
        n_snps += got; n_missing += miss
        del ages, a0_list, a1_list
        print(f"  [{name}] chr{chrom}: {got} used, {miss} missing", flush=True)
    frac = survival / survival[0]
    np.savetxt(f"{OUT}/singer_age_dist.{name}.txt",
               np.column_stack([GRID, GRID * 29 / 1e6, survival, frac * 100]),
               header="age_gen\tage_mya\texpected_snps\tpct_older", delimiter="\t", comments="")
    print(f"[{name}] N_snps={n_snps:,}  missing={n_missing:,}  ({maf_note})", flush=True)
    return n_snps

if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    print("Loading SNP sets...", flush=True)
    shared = load_shared_positions()
    common = load_common_positions(maf_min=0.05)
    ns = sum(len(v) for v in shared.values()); nc = sum(len(v) for v in common.values())
    print(f"shared positions: {ns:,} | common (MAF>=5%) positions: {nc:,}", flush=True)
    run_set("shared", shared, "SINGER-dated shared SNPs")
    run_set("all_common", common, "1000GP MAF>=5% with SINGER age")
    print("DONE")
