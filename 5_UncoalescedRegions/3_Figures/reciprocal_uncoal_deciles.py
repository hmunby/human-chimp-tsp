#!/usr/bin/env python3
"""
Two-panel figure:
but:
  Panel A: % of segments *reciprocally* uncoalesced at 5.5 Mya in BOTH humans AND chimpanzees,
           by BMAP bin (bin 1 = strongest background selection ... bin N = weakest).
  Panel B: number of *new* TSPs per BMAP bin.

Produced at two resolutions: BMAP deciles (N=10) and BMAP quartiles (N=4).

Reciprocal-uncoalesced fraction per bin (whole genome, MHC excluded):
   numerator   = length of the human-uncoal x chimp-uncoal intersection
                 (expected_uncoal_intersect.whole_genome.bed, overlap>0)
   denominator = length of the human-callable x chimp-callable intersection
                 (all_segments_intersect.whole_genome.bed)
   frac(bin) = sum(numerator length in bin) / sum(denominator length in bin)

Both numerator and denominator live in the human n chimp intersect space and are binned by
their (human-side) BMAP value into genome-wide BMAP quantile bins, so the ratio is
"of the sequence jointly callable in both species, what fraction is jointly uncoalesced".
SE is binomial on the per-bin denominator segment count.

BMAP bin edges are genome-wide BMAP quantiles from the 1kb windowed BMAP bed; for N=10 these
reproduce the canonical decile cutoffs. Reproducible: re-run pointing --tsps at the final set.
"""
import argparse
import os
import subprocess
import tempfile
import os
from collections import Counter
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

APRIL = "resources"
BMAP  = "resources/bmap"
NUM   = f"{APRIL}/expected_uncoal/intersect/expected_uncoal_intersect.whole_genome.bed"  # -wao full 
DEN   = f"{APRIL}/all_segments/intersect/all_segments_intersect.whole_genome.bed"           # chrom,start,end,phylop,bmap,len
DECILE_BED = f"{BMAP}/bmap_files/windowed/whole_genome.1kb.bmap_deciles.bed"
CHAIN = "resources/b37tohg38.nochr.over.chain"
LIFTOVER = "liftOver"
OUTDIR = "results"
# Per-variant TSP gene assignments (carries BMAP + assigned_gene), used by the gene panel.
# --tsps overrides this too, so every panel reflects the same TSP set.
VARGENE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "../../4_SINGER_GenealogicalReconstruction_TSP/4_TSPRegions/"
                       "3_GeneAssignment/variant_gene_assignments.pooled.tsv")
# MHC (hg38), excluded throughout
MHC_CHR, MHC_S, MHC_E = "chr6", 28510120, 33480577

# Project convention: all TSP data is drawn in standard bright red.
TSP_BAR  = "#FF0000"   # bright red bar fill
TSP_EDGE = "#FF0000"   # bar edge (single TSP red throughout, matches Figs 2/3)
TSP_DARK = "#FF0000"   # count labels
TSP_LINE = "#FF0000"   # density line

_BMAP_CACHE = None


def bmap_bin_edges(nbins):
    """Interior BMAP quantile edges (len nbins-1) from the genome-wide 1kb windowed BMAP."""
    global _BMAP_CACHE
    if _BMAP_CACHE is None:
        vals = []
        for line in open(DECILE_BED):
            f = line.split("\t")
            try:
                vals.append(float(f[3]))
            except ValueError:
                pass
        _BMAP_CACHE = np.array(vals)
    q = np.quantile(_BMAP_CACHE, np.linspace(0, 1, nbins + 1))
    return q[1:-1]


def to_bin(b, edges):
    # bin 1 = b < edges[0] ; ... ; bin N = b >= edges[-1]
    for i, u in enumerate(edges):
        if b < u:
            return i + 1
    return len(edges) + 1


def in_mhc(chrom, start, end):
    return chrom == MHC_CHR and start < MHC_E and end > MHC_S


def reciprocal_fractions(edges):
    nbins = len(edges) + 1
    num_len = {d: 0.0 for d in range(1, nbins + 1)}
    den_len = {d: 0.0 for d in range(1, nbins + 1)}
    den_n   = {d: 0 for d in range(1, nbins + 1)}

    # numerator: expected_uncoal_intersect.whole_genome.bed (13 cols, -wao)
    #   0 chrom_h 1 start_h 2 end_h 3 p0_h 4 phylop_h 5 bmap_h ... 12 overlap_len
    for line in open(NUM):
        p = line.rstrip("\n").split("\t")
        ov = int(p[12])
        if ov <= 0:
            continue
        if p[5] == "NA":
            continue
        if in_mhc(p[0], int(p[1]), int(p[2])):
            continue
        num_len[to_bin(float(p[5]), edges)] += ov

    # denominator: all_segments_intersect.whole_genome.bed (6 cols)
    #   0 chrom 1 start 2 end 3 phylop 4 bmap 5 length
    for line in open(DEN):
        p = line.rstrip("\n").split("\t")
        if p[4] == "NA":
            continue
        if in_mhc(p[0], int(p[1]), int(p[2])):
            continue
        d = to_bin(float(p[4]), edges)
        den_len[d] += int(p[5])
        den_n[d] += 1

    rows = []
    for d in range(1, nbins + 1):
        frac = num_len[d] / den_len[d] if den_len[d] else 0.0
        n = max(den_n[d], 1)
        se = np.sqrt(frac * (1 - frac) / n)
        rows.append((d, frac, se, num_len[d], den_len[d], den_n[d]))
    return pd.DataFrame(rows, columns=["bin", "frac", "se", "num_len", "den_len", "den_n"])


def _load_records(edges):
    """Per-record (bin, value, chrom-code) arrays for the numerator (overlap bp) and
    denominator (segment bp), MHC excluded. Used to bootstrap the Panel A fraction."""
    ccode = {}
    def code(c):
        return ccode.setdefault(c, len(ccode))
    nb, nv, ncc = [], [], []
    for line in open(NUM):
        p = line.rstrip("\n").split("\t")
        ov = int(p[12])
        if ov <= 0 or p[5] == "NA" or in_mhc(p[0], int(p[1]), int(p[2])):
            continue
        nb.append(to_bin(float(p[5]), edges)); nv.append(ov); ncc.append(code(p[0]))
    db, dv, dcc = [], [], []
    for line in open(DEN):
        p = line.rstrip("\n").split("\t")
        if p[4] == "NA" or in_mhc(p[0], int(p[1]), int(p[2])):
            continue
        db.append(to_bin(float(p[4]), edges)); dv.append(int(p[5])); dcc.append(code(p[0]))
    return (np.array(nb, np.int32), np.array(nv, np.float64), np.array(ncc, np.int32),
            np.array(db, np.int32), np.array(dv, np.float64), np.array(dcc, np.int32),
            len(ccode))


def bootstrap_ci(edges, method="segment", B=200, seed=0):
    """Bootstrap 2.5/97.5% CI for the per-bin reciprocal-uncoalesced fraction.
      method='segment'    -> resample individual segments with replacement (matches
                             Figure 6; numerator and denominator resampled independently).
      method='chromosome' -> resample whole chromosomes with replacement (block bootstrap;
                             accounts for spatial autocorrelation / linkage between segments).
    Returns (lo, hi) per bin, as fractions."""
    nbins = len(edges) + 1
    nb, nv, ncc, db, dv, dcc, C = _load_records(edges)
    rng = np.random.default_rng(seed)

    def bins(bidx, val):
        return np.bincount(bidx, weights=val, minlength=nbins + 1)[1:nbins + 1]

    fr = np.empty((B, nbins))
    if method == "segment":
        Nn, Nd = len(nv), len(dv)
        for b in range(B):
            ni = rng.integers(0, Nn, Nn); di = rng.integers(0, Nd, Nd)
            num = bins(nb[ni], nv[ni]); den = bins(db[di], dv[di])
            fr[b] = np.where(den > 0, num / den, np.nan)
    elif method == "chromosome":
        num_cb = np.zeros((C, nbins)); den_cb = np.zeros((C, nbins))
        np.add.at(num_cb, (ncc, nb - 1), nv)
        np.add.at(den_cb, (dcc, db - 1), dv)
        for b in range(B):
            pick = rng.integers(0, C, C)
            num = num_cb[pick].sum(axis=0); den = den_cb[pick].sum(axis=0)
            fr[b] = np.where(den > 0, num / den, np.nan)
    else:
        raise ValueError(method)
    return np.nanpercentile(fr, 2.5, axis=0), np.nanpercentile(fr, 97.5, axis=0)


def tsp_bins(tsp_file, edges):
    """Lift (hg19/b37) TSP positions to hg38 and count per BMAP bin."""
    nbins = len(edges) + 1
    df = pd.read_csv(tsp_file, sep="\t")
    if "CHROM" not in df.columns:
        raise SystemExit(f"{tsp_file} has no CHROM/POS columns. The count panel needs an "
                         "hg19/b37 TSP table; the gene panel's assignment table (which is "
                         "hg38) goes to --vargene instead.")
    counts = {d: 0 for d in range(1, nbins + 1)}
    with tempfile.TemporaryDirectory() as td:
        b37 = os.path.join(td, "tsp.b37.bed")
        h38 = os.path.join(td, "tsp.hg38.bed")
        un  = os.path.join(td, "tsp.un.bed")
        with open(b37, "w") as f:
            for _, r in df.iterrows():
                c = str(r["CHROM"]); pos = int(r["POS"])
                f.write(f"{c}\t{pos-1}\t{pos}\t{c}_{pos}\n")
        subprocess.run([LIFTOVER, b37, CHAIN, h38, un], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        n_unmapped = (sum(1 for _ in open(un)) // 2) if os.path.exists(un) else 0
        lifted = os.path.join(td, "tsp.hg38.chr.bed")
        with open(h38) as fi, open(lifted, "w") as fo:
            for line in fi:
                c, s, e, name = line.split("\t")[:4]
                c = c if c.startswith("chr") else "chr" + c
                fo.write(f"{c}\t{s}\t{e}\t{name}")
        out = subprocess.run(
            ["bedtools", "intersect", "-a", lifted, "-b", DECILE_BED, "-wa", "-wb"],
            capture_output=True, text=True, check=True).stdout
    seen = set()
    for line in out.strip().split("\n"):
        if not line:
            continue
        f = line.split("\t")
        name = f[3]
        bmap = float(f[-2])   # windowed bed: ... chrom start end bmap decile
        if name in seen:
            continue
        seen.add(name)
        counts[to_bin(bmap, edges)] += 1
    return counts, len(df), len(seen), n_unmapped


def genes_by_bin(edges):
    """Per BMAP bin, list of (gene, n_TSPs) from the per-variant TSP gene assignments,
    sorted by descending SNP count then gene name."""
    nbins = len(edges) + 1
    d = pd.read_csv(VARGENE, sep="\t")
    # restrict to the candidate (860-anchored) TSP set, matching tsp_snps_final.txt
    keyf = f"{OUTDIR}/results/candidate_tsp_hg38_keys.txt"
    if os.path.exists(keyf):
        keys = {l.strip() for l in open(keyf) if l.strip()}
        d = d[[f"{c}_{p}" in keys for c, p in
               zip(d["CHROM_HG38"].astype(str), d["POSITION_HG38"].astype(int))]]
    d["BMAP"] = pd.to_numeric(d["BMAP"], errors="coerce")
    out = {b: [] for b in range(1, nbins + 1)}
    for b in range(1, nbins + 1):
        lo = -1e18 if b == 1 else edges[b - 2]
        hi = 1e18 if b == nbins else edges[b - 1]
        sub = d[(d["BMAP"] >= lo) & (d["BMAP"] < hi)]
        vc = sub["assigned_gene"].dropna().astype(str).value_counts()
        out[b] = sorted(vc.items(), key=lambda kv: (-kv[1], kv[0]))
    return out


def make_figure(nbins, label, frac_df, counts, out, mode="density", genes=None):
    if genes is None:
        fig = plt.figure(figsize=(12, 7), dpi=300, constrained_layout=True)
        spec = gridspec.GridSpec(ncols=1, nrows=2, figure=fig, height_ratios=[3, 1])
        axc = None
    else:
        # gene-panel height auto-scales with the max number of genes in any bin, so the
        # per-gene row spacing stays constant/legible. Panel heights are set in inches
        # (height_ratios == inch values, so the figure size preserves them exactly).
        maxg = max((len(g) for g in genes.values()), default=1)
        A_in, B_in = 2.8, 1.8                        # taller A so the 2-line y-label fits; taller B
        C_in = max(maxg * 0.26, 1.2)                 # ~0.26 in per gene row
        figH = A_in + B_in + C_in
        fig = plt.figure(figsize=(12, figH), dpi=300, constrained_layout=True)
        fig.set_constrained_layout_pads(hspace=0.08, h_pad=0.1)   # separate stacked panels/labels
        spec = gridspec.GridSpec(ncols=1, nrows=3, figure=fig,
                                 height_ratios=[A_in, B_in, C_in])
        axc = fig.add_subplot(spec[2, 0])
    ax  = fig.add_subplot(spec[0, 0])
    ax2 = fig.add_subplot(spec[1, 0])

    perc = frac_df["frac"] * 100
    perc_se = frac_df["se"] * 100
    # Error bars: bootstrap 2.5/97.5% CI when provided (matches Figure 6's method),
    # else fall back to the binomial SE. Marker/line weights match Figure 6.
    if "ci_lo" in frac_df.columns:
        pv = perc.to_numpy()
        lo = frac_df["ci_lo"].to_numpy() * 100; hi = frac_df["ci_hi"].to_numpy() * 100
        yerr = np.clip(np.vstack([pv - lo, hi - pv]), 0, None); ytop = float(hi.max())
    else:
        yerr = perc_se; ytop = float((perc + perc_se).max())
    ax.errorbar(frac_df["bin"], perc, yerr=yerr, fmt='-o', color="#4C72B0",
                markersize=4, capsize=3, elinewidth=1.2)
    ax.text(-0.07, 1.05, 'A', transform=ax.transAxes, fontsize=16, fontweight='bold',
            va='top', ha='right')
    ax.set_ylabel("Percentage of segments\nuncoalesced in both species", fontsize=11, fontweight='bold', labelpad=6)
    # shared x-axis: tick numbers only on the bottom panel (B); hide them here
    ax.set_xticks(range(1, nbins + 1)); ax.tick_params(labelbottom=False)
    ax.set_xlim(0.5, nbins + 0.5)
    ax.set_ylim(0, ytop * 1.15)

    ax2.text(-0.07, 1.15, 'B', transform=ax2.transAxes, fontsize=16, fontweight='bold',
             va='top', ha='right')
    bins = list(range(1, nbins + 1))
    # raw counts as bars (left axis) -- all TSP data in red
    cnt = [counts[d] for d in bins]
    ax2.bar(bins, cnt, width=0.65, color=TSP_BAR, edgecolor=TSP_EDGE, zorder=2)
    for d, c in zip(bins, cnt):
        ax2.text(d, c + max(cnt) * 0.02, str(c), ha='center', va='bottom',
                 color=TSP_DARK, fontsize=9, fontweight='bold', zorder=3)
    ax2.set_ylabel("Number of TSPs", color="black", labelpad=8, fontweight='bold')
    ax2.tick_params(axis='y', colors="black")
    ax2.set_ylim(0, max(cnt) * 1.28 + 1)
    ax2.set_xticks(bins); ax2.set_xlim(0.5, nbins + 0.5)
    xtl = [str(b) for b in bins]
    xtl[0] = f"{bins[0]}\n" + r"$\mathbf{(strongest)}$"
    xtl[-1] = f"{bins[-1]}\n" + r"$\mathbf{(weakest)}$"
    ax2.set_xticklabels(xtl)
    ax2.set_xlabel(f"Background selection {label}", fontsize=12, fontweight='bold')
    if mode == "density":
        ax2.set_zorder(1); ax2.patch.set_visible(False)
        # density (TSPs per Mb of reciprocally-uncoalesced sequence) on twin axis -- TSP data, red
        dens = [counts[d] / (fr.num_len / 1e6) if fr.num_len > 0 else np.nan
                for d, fr in zip(bins, frac_df.itertuples(index=False))]
        ax2b = ax2.twinx()
        ax2b.plot(bins, dens, '-o', color=TSP_LINE, markersize=6, lw=1.5, zorder=4)
        ax2b.set_ylabel("TSPs per Mb reciprocally-\nuncoalesced sequence", color=TSP_LINE,
                        fontsize=9, labelpad=8)
        ax2b.tick_params(axis='y', colors=TSP_LINE)
        ax2b.set_ylim(0, np.nanmax(dens) * 1.2)

    # ---- Panel C: assigned gene names per bin (stacked under each column) ----
    if axc is not None:
        axc.text(-0.07, 1.02, 'C', transform=axc.transAxes, fontsize=16, fontweight='bold',
                 va='top', ha='right')
        axc.set_xlim(0.5, nbins + 0.5); axc.set_ylim(0, 1)
        axc.set_xticks([]); axc.set_yticks([])
        for spine in axc.spines.values():
            spine.set_visible(False)
        axc.set_ylabel("Assigned genes", labelpad=8, fontweight='bold')
        maxg = max((len(g) for g in genes.values()), default=1)
        # genes appearing in >1 BMAP bin (SNPs span multiple deciles/quartiles) -> asterisk
        bin_span = Counter(g for b in genes for g, _ in genes[b])
        multi = {g for g, c in bin_span.items() if c > 1}
        for b in range(1, nbins + 1):
            gl = genes.get(b, [])
            for j, (g, n) in enumerate(gl):
                # pretty-print unnamed Ensembl IDs; '*' if multi-bin; count only when >1 SNP
                name = g if not g.startswith("ENSG") else g[:4] + "…" + g[-4:]
                lbl = name + ("*" if g in multi else "") + (f" ({n})" if n > 1 else "")
                axc.text(b, 0.985 - j * (0.97 / max(maxg, 1)), lbl, ha='center', va='top',
                         fontsize=10, fontstyle='italic', color="#333333")
            if not gl:
                axc.text(b, 0.5, "–", ha='center', va='center', fontsize=8, color="#999999")

    # align the left-hand y-axis labels across panels A, B (and C if present)
    _axes = [ax, ax2] + ([axc] if axc is not None else [])
    fig.align_ylabels(_axes)
    fig.savefig(out, dpi=300, bbox_inches='tight')
    print(f"wrote {out}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vargene", help="per-variant gene assignments (hg38) for the gene panel")
    ap.add_argument("--outdir", help="write results/figures here instead of the default")
    ap.add_argument("--tsps", default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "../../6_TSP_Analyses/snp_set_comparisons/misc_files/tsp_snps_pooled.hg19.txt"))
    a = ap.parse_args()
    if a.tsps:
        global VARGENE
    if a.vargene:
        global VARGENE
        VARGENE = a.vargene
    if a.outdir:
        global OUTDIR
        OUTDIR = a.outdir
        import os as _os
    for _sub in ('results', 'figures'):
        _os.makedirs(f'{OUTDIR}/{_sub}', exist_ok=True)

    for nbins, label in [(10, "decile"), (4, "quartile")]:
        edges = bmap_bin_edges(nbins)
        frac_df = reciprocal_fractions(edges)
        counts, n_total, n_assigned, n_unmapped = tsp_bins(a.tsps, edges)
        tag = f"{label}s"
        frac_df.to_csv(f"{OUTDIR}/results/reciprocal_uncoal.{tag}.tsv", sep="\t", index=False)
        dens = [counts[d] / (nl / 1e6) if nl > 0 else float("nan")
                for d, nl in zip(frac_df["bin"], frac_df["num_len"])]
        pd.DataFrame({"bin": list(frac_df["bin"]), "n_tsps": [counts[d] for d in frac_df["bin"]],
                      "tsps_per_Mb_reciprocal_uncoal": dens}).to_csv(
            f"{OUTDIR}/results/tsp_counts.{tag}.tsv", sep="\t", index=False)
        print(f"\n===== BMAP {tag} (edges={[round(x,1) for x in edges]}) =====")
        print(frac_df.to_string(index=False))
        print(f"TSPs: total={n_total}, assigned={n_assigned}, liftover-unmapped={n_unmapped}")
        print("per-bin TSP counts:", counts)
        gbb = genes_by_bin(edges)
        pd.DataFrame([(b, ", ".join(f"{g} ({n})" for g, n in gbb[b])) for b in gbb],
                     columns=["bin", "assigned_genes"]).to_csv(
            f"{OUTDIR}/results/genes_by_bin.{tag}.tsv", sep="\t", index=False)
        # segment bootstrap -> main figures; chromosome-block -> separate ".chromblock" figures
        for method, suffix in [("segment", ""), ("chromosome", ".chromblock")]:
            lo, hi = bootstrap_ci(edges, method=method, B=200, seed=0)
            fdf = frac_df.copy(); fdf["ci_lo"] = lo; fdf["ci_hi"] = hi
            print(f"  {method} bootstrap 95% CI (%): "
                  + ", ".join(f"{int(b)}:[{l*100:.3f},{h*100:.3f}]"
                              for b, l, h in zip(fdf["bin"], lo, hi)))
            make_figure(nbins, label, fdf, counts,
                        f"{OUTDIR}/figures/reciprocal_uncoal_bmap_{tag}.counts{suffix}.png", mode="counts")
            make_figure(nbins, label, fdf, counts,
                        f"{OUTDIR}/figures/reciprocal_uncoal_bmap_{tag}.counts_density{suffix}.png", mode="density")
            make_figure(nbins, label, fdf, counts,
                        f"{OUTDIR}/figures/reciprocal_uncoal_bmap_{tag}.genes{suffix}.png",
                        mode="counts", genes=gbb)


if __name__ == "__main__":
    main()
