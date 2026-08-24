#!/usr/bin/env python3
"""Build the human expected-uncoalesced segments from the bundled Cousins et al. decode.

    per-individual P(uncoalesced)  ->  sum over the 7 individuals > 1  ->  drop segmental
    duplications  ->  mean phyloP  ->  BMAP  ->  the "unconstrained" cutoffs

This is the human counterpart of build_chimp_and_intersect.py, and writes the human half of
`expected_uncoal/` that that script and the Figure 4/6 scripts read. Both halves must exist
before phylop_annotate_and_gene_test.sh will run.

The 5.5 Mya column. The decode gives P(not yet coalesced) at each of 19 time boundaries; the
one used throughout is PO_5589889, the first boundary past the 5.5 Mya human-chimpanzee split
on the human grid (mu = 1.25e-8, generation time 29). The chimpanzee side uses its own grid and
lands on 5,755,625 years instead -- see ../../1_ChimpPSMC/4_MinimalDecode/.

Segmental duplications are removed from the human segments only, matching the published
analysis; the chimpanzee segments are not SD-filtered.

    python build_human_uncoalesced.py

Writes into ../expected_uncoal/human/ (per chromosome and genome-wide).
"""
import argparse
import gzip
import os
import subprocess
import tempfile

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))

# ---- inputs (edit) -------------------------------------------------------------------------
DATA = os.path.join(HERE, "data")                       # the bundled per-individual columns
OUT = os.path.join(HERE, "..", "expected_uncoal", "human")

# T2T segmental duplications lifted to hg38, bare chromosome names. Not bundled -- see RESOURCES.
SD_BED = os.path.join(HERE, "../resources/chm13v2.0_SD.pos_only.hg38_coords_only.bed")
# UCSC hg38 phyloP 100-way bigWig (~9.9 GB) and the genome-wide 1 kb BMAP decile bed.
BIGWIG = os.path.join(HERE, "../resources/hg38.phyloP100way.bw")
BMAP_DECILE = os.path.join(HERE, "../resources/whole_genome.1kb.bmap_deciles.bed")

BEDTOOLS = "bedtools"
BIGWIGAVG = "bigWigAverageOverBed"
# --------------------------------------------------------------------------------------------

AGE = 5589889          # the PO_ boundary this column represents, years
SCALE = 1_000_000      # the bundled integers are round(P * SCALE)
THRESHOLD = 1.0        # keep windows whose summed probability exceeds this
CUTOFFS = [0.090083, 0.150087, 0.331998, 0.617324]      # 70/80/90/95th phyloP percentiles

# genome-wide concat order: string sort of chr1..chr22, as the pipeline built it
CHR_ORDER = sorted(range(1, 23), key=str)


def load(sample):
    """One individual's column: {chrom: int array}, positions implicit at 1 kb from 0."""
    counts, values = [], []
    with gzip.open(os.path.join(DATA, f"{sample}.PO_5589889.txt.gz"), "rt") as fh:
        for line in fh:
            if line.startswith("#chrom"):
                _, chrom, n = line.split()
                counts.append((int(chrom), int(n)))
            elif not line.startswith("#"):
                values.append(int(line))
    arr = np.array(values, dtype=np.int64)
    out, i = {}, 0
    for chrom, n in counts:
        out[chrom] = arr[i:i + n]
        i += n
    assert i == len(arr), "window counts in the header do not match the number of values"
    return out


def summed_segments(samples):
    """Windows where the summed probability exceeds THRESHOLD, per chromosome.

    Summation is in the scaled integers, so the comparison is exact rather than dependent on
    floating-point accumulation order.
    """
    cols = {s: load(s) for s in samples}
    cut = int(round(THRESHOLD * SCALE))
    for chrom in range(1, 23):
        n = max(len(cols[s][chrom]) for s in samples)
        total = np.zeros(n, dtype=np.int64)
        for s in samples:
            v = cols[s][chrom]
            total[:len(v)] += v
        idx = np.flatnonzero(total > cut)
        yield chrom, idx * 1000, total[idx] / SCALE


def bmap_lookup():
    """(chrom, rounded start) -> raw B value."""
    bd = pd.read_csv(BMAP_DECILE, sep="\t", header=None, names=["c", "s", "e", "bmap", "dec"])
    return dict(zip(zip(bd["c"], bd["s"]), bd["bmap"]))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--skip-phylop", action="store_true",
                    help="stop after the SD-filtered segments (no bigWig needed)")
    a = ap.parse_args()

    samples = sorted(f.split(".")[0] for f in os.listdir(DATA) if f.endswith(".txt.gz"))
    print(f"  individuals: {len(samples)} ({', '.join(samples)})")

    phydir = os.path.join(a.out, "phylop")
    os.makedirs(phydir, exist_ok=True)
    bmap = None if a.skip_phylop else bmap_lookup()
    kept_total = 0
    per_chrom = {}

    with tempfile.TemporaryDirectory() as tmp:
        for chrom, starts, probs in summed_segments(samples):
            # 1. summed windows, bare chromosome names to match the SD mask
            raw = os.path.join(tmp, f"c{chrom}.bed")
            with open(raw, "w") as fh:
                fh.writelines(f"{chrom}\t{s}\t{s + 1000}\t{p!r}\n" for s, p in zip(starts, probs))

            # 2. drop any window touching a segmental duplication (-v drops on any overlap)
            sd_out = os.path.join(a.out, f"expected_uncoal_human_{AGE}.chr{chrom}.bed")
            with open(sd_out, "w") as fh:
                subprocess.run([BEDTOOLS, "intersect", "-a", raw, "-b", SD_BED, "-v"],
                               stdout=fh, check=True)
            n = sum(1 for _ in open(sd_out))
            kept_total += n
            print(f"  chr{chrom:<2} {len(starts):>6,} windows > {THRESHOLD} "
                  f"-> {n:>6,} after SD removal", flush=True)
            if a.skip_phylop:
                continue

            # 3. mean phyloP per window. The name column is the row number, not the probability:
            #    bigWigAverageOverBed rejects duplicate names, and rounded probabilities do
            #    collide. -bedOut preserves row order, so the probability pastes back after.
            pre = os.path.join(tmp, f"c{chrom}.chr.bed")
            subprocess.run(f"awk 'BEGIN{{OFS=\"\\t\"}} {{print \"chr\"$1,$2,$3,NR}}' {sd_out} > {pre}",
                           shell=True, check=True)
            scored = os.path.join(tmp, f"c{chrom}.phylop.bed")
            subprocess.run([BIGWIGAVG, BIGWIG, pre, os.path.join(tmp, f"c{chrom}.tab"),
                            f"-bedOut={scored}"], check=True, capture_output=True)

            # 4. probability back into the name column, then BMAP from the 1 kb decile bed
            phy = os.path.join(phydir, f"expected_uncoal_human_{AGE}.chr{chrom}.phylop100way.bed")
            out6 = os.path.join(phydir,
                                f"expected_uncoal_human_{AGE}.chr{chrom}.phylop100way.bmap.bed")
            probs = [l.split("\t")[3].rstrip("\n") for l in open(sd_out)]
            with open(scored) as fin, open(phy, "w") as f5, open(out6, "w") as f6:
                for i, line in enumerate(fin):
                    c, s, e, _, ph = line.rstrip("\n").split("\t")
                    f5.write(f"{c}\t{s}\t{e}\t{probs[i]}\t{ph}\n")
                    f6.write(f"{c}\t{s}\t{e}\t{probs[i]}\t{ph}\t"
                             f"{bmap.get((c, 1000 * round(int(s) / 1000)), 'NA')}\n")
            per_chrom[chrom] = out6

    print(f"\n  genome-wide: {kept_total:,} windows")
    if a.skip_phylop:
        return

    # 5. genome-wide concat, then the four "unconstrained" phyloP cutoffs
    gw = os.path.join(phydir, f"expected_uncoal_human_{AGE}.genome_wide.phylop100way.bmap.bed")
    with open(gw, "w") as out:
        for chrom in CHR_ORDER:
            with open(per_chrom[chrom]) as fh:
                out.writelines(fh)
    print(f"  wrote {gw} ({sum(1 for _ in open(gw)):,} rows)")

    for c in CUTOFFS:
        f = os.path.join(phydir,
                         f"expected_uncoal_human_{AGE}.genome_wide.phylop100way."
                         f"unconstrained_cutoff_{c:.6f}.bed")
        with open(f, "w") as out:
            subprocess.run(["awk", "-F", "\t", f"$5 < {c}", gw], stdout=out, check=True)
        print(f"  cutoff {c:.6f}: {sum(1 for _ in open(f)):,} rows")


if __name__ == "__main__":
    main()
