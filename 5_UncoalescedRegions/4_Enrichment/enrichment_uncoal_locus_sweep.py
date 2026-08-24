#!/usr/bin/env python3
"""
Merged-region enrichment for uncoalesced regions, swept over merge distance D.

For each D in the sweep:
  * merge regions whose gap <= D into loci; span = min(start)..max(end) including internal gaps
  * candidate genes of a locus = UNION of its regions' closest-gene sets
  * locus BMAP = length-weighted mean over constituent regions -> quartile, as merge_bed_weighted
  * ONE gene per locus, drawn at random, averaged over N_TIEBREAK resolutions
  * null draws one SPAN-MATCHED contiguous background block per locus, from the same quartile.
  * Fraction statistic for fold and permutation p, BH within each tie-break resolution.

Usage: python3 enrichment_uncoal_locus_sweep.py [NPERM] [N_TIEBREAK] [D1,D2,...]
       defaults: 10000 200 0,5000,10000,25000
"""

from pathlib import Path
HERE = Path(__file__).resolve().parent
# Repo-relative paths.
AP = str(HERE / "resources")
GA = str(HERE / "../../6_TSP_Analyses/enrichment/geneset")
GL = str(HERE / "../../6_TSP_Analyses/enrichment/glycoprotein")
import glob, sys
import numpy as np
import pandas as pd
from scipy.stats import hypergeom




MHC = ("chr6", 28510120, 33480577)

NPERM = int(sys.argv[1]) if len(sys.argv) > 1 else 10000
NTIE = int(sys.argv[2]) if len(sys.argv) > 2 else 200
DLIST = [int(x) for x in sys.argv[3].split(",")] if len(sys.argv) > 3 else [0, 5000, 10000, 25000]
# Which D values get the (expensive) null. Others report observed-side locus/gene counts only.
NULL_FOR = set(int(x) for x in sys.argv[4].split(",")) if len(sys.argv) > 4 else set(DLIST)
# Optional cap on region->gene distance (bp). A region whose closest gene is further than this is
# DROPPED entirely, and so is any null block whose closest gene is further than this, so the two
# sides stay comparable. Default: no cap.
# Motivation: 4 of the 28 Q1 membrane-glycoprotein genes sit 508kb-1.08Mb away (NRG3, UNC5D,
# TMEM182, CADM2), all large genes in gene-poor neighbourhoods where "nearest gene" is weak.
DISTCAP = float(sys.argv[5]) if len(sys.argv) > 5 else np.inf
MINKS = [3, 2]
GLYCO = ["KW-0325 glycoprotein", "membrane glycoprotein (KW-0325 & KW-0472)"]


def load_gmt(p):
    d = {}
    for line in open(p):
        f = line.rstrip("\n").split("\t")
        if len(f) > 2:
            d[f[0]] = {g.upper() for g in f[2:] if g}
    return d


def load_tsv(p):
    gl = pd.read_csv(p, sep="\t")
    c = [x for x in gl.columns if x.lower() in
         ("hgnc_symbol", "gene_name_primary", "symbol", "gene names (primary)")][0]
    return {str(x).upper() for x in gl[c].dropna()}


def bh(p):
    p = np.asarray(p, float); m = len(p)
    if m == 0: return p
    o = np.argsort(p)
    r = np.minimum.accumulate((p[o] * m / (np.arange(m) + 1))[::-1])[::-1]
    out = np.empty(m); out[o] = np.clip(r, 0, 1); return out


CATS = {}
for g in sorted(glob.glob(f"{GA}/genesets/*.gmt")): CATS.update(load_gmt(g))
CATS[GLYCO[0]] = load_tsv(f"{GL}/human_glycoproteins_KW0325_with_coords.tsv")
CATS[GLYCO[1]] = load_tsv(f"{GL}/human_membrane_glycoproteins_KW0325_KW0472_with_coords.tsv")

bc = pd.read_csv(f"{AP}/all_segments/intersect/all_segments_intersect.whole_genome.bmap.length_counts.txt",
                 header=None, sep=" ", names=["bmap", "count"]).sort_values("bmap")
bc = bc[bc.bmap != "NA"].copy(); bc["bmap"] = bc.bmap.astype(float)
bc["cf"] = bc["count"].cumsum() / bc["count"].sum()
EDGES = np.array([bc[bc.cf >= q].bmap.iloc[0] for q in [.25, .5, .75]])
qv = lambda b: np.where(np.isnan(np.asarray(b, float)), 0,
                        np.searchsorted(EDGES, np.asarray(b, float), side="left") + 1).astype(np.int8)

# ------------------------------------------------------------------ observed regions
C = ["chrom", "start", "end", "p0", "phylop", "bmap", "cg", "gs", "ge", "gene", "distance"]
d = pd.read_csv(f"{AP}/expected_uncoal/intersect/expected_uncoal_intersect.whole_genome.exact.merged.closest_genes.bed",
                header=None, sep="\t", names=C)
d["bmap"] = pd.to_numeric(d.bmap, errors="coerce"); d = d.dropna(subset=["bmap"])
d = d[~((d.chrom == MHC[0]) & (d.start < MHC[2]) & (d.end > MHC[1]))]
d["gene"] = d.gene.astype(str).str.upper()
if np.isfinite(DISTCAP):
    d["distance"] = pd.to_numeric(d.distance, errors="coerce")
    n0 = d.groupby(["chrom", "start", "end"]).ngroups
    d = d[d.distance <= DISTCAP]
    print(f"DIST CAP {DISTCAP:,.0f} bp: {n0} regions -> "
          f"{d.groupby(['chrom','start','end']).ngroups} with a gene within the cap", flush=True)
reg = d.groupby(["chrom", "start", "end"], as_index=False).agg(bmap=("bmap", "first"))
reg = reg.sort_values(["chrom", "start"]).reset_index(drop=True)
cand = d.groupby(["chrom", "start", "end"]).gene.apply(lambda s: sorted(set(s))).to_dict()
print(f"{len(reg)} bookended-merged regions", flush=True)

# ------------------------------------------------------------------ background
BC = ["chrom", "start", "end", "phylop", "bmap", "ov", "cg", "gs", "ge", "gene", "distance"]
bg = pd.read_csv(f"{AP}/all_segments/intersect/all_segments_intersect.whole_genome.closest_genes.dedup.bed",
                 header=None, sep="\t", names=BC,
                 usecols=["chrom", "start", "end", "bmap", "gene", "distance"])
bg = bg.sort_values(["chrom", "start"]).reset_index(drop=True)
bg["bmap"] = pd.to_numeric(bg.bmap, errors="coerce")
bg["distance"] = pd.to_numeric(bg.distance, errors="coerce").fillna(1 << 30)
bg["gene"] = bg.gene.astype(str).str.upper()
inmhc = ((bg.chrom == MHC[0]) & (bg.start < MHC[2]) & (bg.end > MHC[1])).values
NB = len(bg)
print(f"{NB} background segments", flush=True)

VOCAB = pd.Index(sorted(set(bg.gene) | set(d.gene)))
NG = len(VOCAB)
bgcode = VOCAB.get_indexer(bg.gene.values)
CIDX = {c: VOCAB.get_indexer(sorted(g)) for c, g in CATS.items()}
CIDX = {c: i[i >= 0] for c, i in CIDX.items()}

ccode, cuniq = pd.factorize(bg.chrom.values)
OFF = (np.arange(len(cuniq), dtype=np.int64) + 1) * 10**10       # monotonic global coords
gstart = OFF[ccode] + bg.start.values.astype(np.int64)
gend = OFF[ccode] + bg.end.values.astype(np.int64)
dist = bg.distance.values.astype(np.float64)
blen = (bg.end.values - bg.start.values).astype(np.float64)
bbm = bg.bmap.values.astype(np.float64)

adj = np.zeros(NB, dtype=bool)
adj[:-1] = (ccode[1:] == ccode[:-1]) & (bg.start.values[1:] == bg.end.values[:-1]) & ~inmhc[1:] & ~inmhc[:-1]
newrun = np.ones(NB, dtype=bool); newrun[1:] = ~adj[:-1]
grp = np.cumsum(newrun) - 1
rs = np.flatnonzero(newrun); re_ = np.concatenate([rs[1:] - 1, [NB - 1]])
RUNEND = re_[grp]                                                 # last index of this segment's run
ok = ~np.isnan(bbm)
W = np.concatenate([[0.0], np.cumsum(np.where(ok, bbm * blen, 0.0))])
DD = np.concatenate([[0.0], np.cumsum(np.where(ok, blen, 0.0))])


def _closest(s, lens):
    """Closest gene and its distance for blocks starting at s of length lens.
    Bucketed padded gather: rectangles per length band so short blocks don't pay for the longest."""
    gene = np.empty(len(s), dtype=np.int64); dmin = np.full(len(s), np.inf)
    order = np.argsort(lens)
    for lo, hi in [(1, 4), (5, 16), (17, 64), (65, 1 << 30)]:
        sel = order[(lens[order] >= lo) & (lens[order] <= hi)]
        if len(sel) == 0: continue
        ml = int(lens[sel].max())
        ii = s[sel][:, None] + np.arange(ml)[None, :]
        m = np.arange(ml)[None, :] < lens[sel][:, None]
        v = np.where(m, dist[np.minimum(ii, NB - 1)], np.inf)
        am = v.argmin(axis=1)
        gene[sel] = bgcode[s[sel] + am]
        dmin[sel] = v[np.arange(len(sel)), am]
    return gene, dmin


def draw_null_genes(spans, q, rng, tries=80):
    """One span-matched, quartile-matched contiguous background block per locus -> its closest gene.

    When DISTCAP is set, a block whose closest gene exceeds the cap is REJECTED AND REDRAWN, not
    discarded. Discarding would leave the null with fewer units than the observed set, which is
    conditioned to have a gene within the cap by construction. (Discarding was the original bug:
    it drove the null/observed gene-set ratio to 0.66-0.78 in Q2/Q3.)
    """
    L = len(spans)
    out = np.full(L, -1, dtype=np.int64)
    todo = np.ones(L, dtype=bool)
    for _ in range(tries):
        idxs = np.flatnonzero(todo)
        if len(idxs) == 0: break
        cs = rng.integers(0, NB, size=len(idxs))
        cj = np.searchsorted(gend, gstart[cs] + spans[idxs], side="left")
        good = (cj < NB) & (cj <= RUNEND[cs])
        cj = np.minimum(cj, NB - 1)
        den = DD[cj + 1] - DD[cs]; num = W[cj + 1] - W[cs]
        with np.errstate(invalid="ignore", divide="ignore"):
            bq = np.where(den > 0, num / np.where(den > 0, den, 1), np.nan)
        good &= (qv(bq) == q)
        gi = np.flatnonzero(good)
        if len(gi) == 0: continue
        g, dm = _closest(cs[gi], (cj[gi] - cs[gi] + 1).astype(np.int64))
        keep = dm <= DISTCAP                       # inf cap -> always true
        sel = idxs[gi[keep]]
        out[sel] = g[keep]
        todo[sel] = False
    return out[out >= 0]


results = {}
for D in DLIST:
    same = reg.chrom.values[1:] == reg.chrom.values[:-1]
    gap = np.where(same, reg.start.values[1:] - reg.end.values[:-1], np.inf)
    newl = np.ones(len(reg), dtype=bool); newl[1:] = (~same) | (gap > D)
    reg["lid"] = np.cumsum(newl) - 1
    rl = (reg.end - reg.start).values.astype(float)
    loci = []
    for lid, g in reg.groupby("lid"):
        w = (g.end - g.start).values.astype(float)
        cs = sorted({x for _, r in g.iterrows() for x in cand[(r.chrom, r.start, r.end)]})
        loci.append((int(qv(np.average(g.bmap.values, weights=w))),
                     int(g.end.max() - g.start.min()), cs))
    LQ = np.array([x[0] for x in loci]); LS = np.array([x[1] for x in loci])
    print(f"\n########## D = {D:,} bp : {len(loci)} loci "
          f"(from {len(reg)} regions) ##########", flush=True)

    for q in [1, 2, 3, 4]:
        sel = np.flatnonzero(LQ == q)
        spans = LS[sel]
        cl = [loci[i][2] for i in sel]
        flat = VOCAB.get_indexer(np.concatenate([np.array(c) for c in cl]))
        cnt = np.array([len(c) for c in cl]); off = np.concatenate([[0], np.cumsum(cnt)[:-1]])

        trng = np.random.default_rng(500 + q)
        obs = np.zeros((NTIE, NG), dtype=bool); nt = np.zeros(NTIE, dtype=np.int32)
        for t in range(NTIE):
            pick = flat[off + trng.integers(0, cnt)]
            u = np.unique(pick); obs[t, u] = True; nt[t] = len(u)

        if D not in NULL_FOR:
            print(f"  Q{q}: {len(sel)} loci, genes {nt.min()}-{nt.max()} (no null run)", flush=True)
            continue

        rng = np.random.default_rng(q)
        masks = np.zeros((NPERM, NG), dtype=bool); sz = np.zeros(NPERM, dtype=np.int32)
        for r in range(NPERM):
            u = np.unique(draw_null_genes(spans, q, rng))
            masks[r, u] = True; sz[r] = len(u)
            if (r + 1) % 5000 == 0: print(f"    Q{q} perm {r+1}/{NPERM}", flush=True)
        print(f"  Q{q}: {len(sel)} loci | obs genes {nt.min()}-{nt.max()} "
              f"| null genes {sz.mean():.0f} (ratio {sz.mean()/nt.mean():.2f})", flush=True)

        # ---- gene-level ("naive") background: every gene that is nearest to ANY callable segment
        # in this BMAP quartile, each counted once. Same nearest-gene procedure as the test side,
        # applied to the whole quartile rather than to the uncoalesced loci. Distance cap applied
        # here too when set, so test and background are filtered the same way.
        # MHC excluded here exactly as on the test side; leaving it in adds ~170 genes to the
        # Q1 background and biases the hypergeometric.
        bsub = bg[(qv(bg.bmap.values) == q) & ~inmhc]
        if np.isfinite(DISTCAP):
            bsub = bsub[bsub.distance <= DISTCAP]
        bgmask = np.zeros(NG, dtype=bool)
        bgmask[VOCAB.get_indexer(bsub.gene.unique())] = True
        nbg = int(bgmask.sum())

        cats, K, P, NF, HF, HP = [], [], [], [], [], []
        for c, idx in CIDX.items():
            if len(idx) == 0: continue
            kt = obs[:, idx].sum(axis=1)
            if kt.max() < 2: continue
            kn = masks[:, idx].sum(axis=1); nf = kn / sz
            ge = NPERM - np.searchsorted(np.sort(nf), kt / nt, side="left")
            kbg = int(bgmask[idx].sum())
            if kbg:
                hf = np.median((kt / nt) / (kbg / nbg))
                hp = float(np.median([hypergeom.sf(a - 1, nbg, kbg, b) for a, b in zip(kt, nt)]))
            else:
                hf, hp = np.inf, 1.0
            cats.append(c); K.append(kt); P.append((1 + ge) / (1 + NPERM)); NF.append(nf.mean())
            HF.append(hf); HP.append(hp)
        K = np.array(K); P = np.array(P); NF = np.array(NF)
        HF = np.array(HF); HP = np.array(HP)
        FOLD = np.median(K / nt[None, :], axis=1) / NF
        print(f"    gene-level background: {nbg} genes in Q{q}", flush=True)
        rows = {}
        for MK in MINKS:
            ps = np.zeros(len(cats)); te = np.zeros(len(cats)); bhv = np.full((len(cats), NTIE), np.nan)
            for t in range(NTIE):
                s2 = np.flatnonzero(K[:, t] >= MK)
                if not len(s2): continue
                a = bh(P[s2, t]); bhv[s2, t] = a; te[s2] += 1; ps[s2] += (a < 0.05)
            with np.errstate(invalid="ignore"):
                rows[MK] = pd.DataFrame(dict(category=cats, k_median=np.median(K, axis=1),
                    fold=FOLD, p_median=np.median(P, axis=1), BH_median=np.nanmedian(bhv, axis=1),
                    frac_tested=te / NTIE, frac_pass=np.where(te > 0, ps / np.maximum(te, 1), 0.0),
                    # gene-level ("naive") hypergeometric against the quartile background
                    fold_genelevel=HF, p_genelevel=HP, BH_genelevel=bh(HP), nbg=nbg))
        results[(D, q)] = rows
        for MK in MINKS:
            r = rows[MK]; rob = r[(r.frac_pass >= .95) & (r.frac_tested >= .95)]
            print(f"    MIN_K={MK}: {len(r)} cats, {len(rob)} robust", flush=True)
            for _, x in r[r.category.isin(GLYCO)].iterrows():
                print(f"       {x.category[:40]:40s} k={x.k_median:5.1f} fold={x.fold:5.2f} "
                      f"BH={x.BH_median:.4f} pass={100*x.frac_pass:3.0f}%", flush=True)

out = []
for (D, q), rows in results.items():
    for MK, r in rows.items():
        r = r.copy(); r["D"] = D; r["quartile"] = q; r["MIN_K"] = MK; out.append(r)
if out:
    # filename carries the parameters, so concurrent/successive runs cannot clobber each other
    tag = "D" + "_".join(str(x // 1000) for x in sorted(NULL_FOR)) + "kb"
    if np.isfinite(DISTCAP):
        tag += f".cap{int(DISTCAP)//1000}kb"
    fn = f"{AP}/results/uncoal_enrichment_locus_sweep.{tag}.tsv"
    pd.concat(out).to_csv(fn, sep="\t", index=False)
    print(f"\nwrote {fn}")
