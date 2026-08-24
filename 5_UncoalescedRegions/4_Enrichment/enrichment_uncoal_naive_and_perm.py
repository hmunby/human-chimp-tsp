"""
Gene-set / glycoprotein enrichment in the reciprocally-uncoalesced regions, per BMAP quartile,
computed TWO ways for EVERY category:

  (1) naive hypergeometric  : gene-level, matched-by-BMAP-quartile background.
        fold_hyperg = (ktest / ntest) / (kbg / nbg)
        p_hyperg    = one-sided hypergeometric P(X >= ktest), population = background genes.
  (2) segment permutation   : resample nseg nearest-gene slots (= observed #uncoalesced segments)
        from the same-quartile all-segment pool, NPERM times.
        fold_perm = ktest / mean(knull);  perm_p = (1 + #{knull >= ktest}) / (1 + NPERM)

Both share the SAME construction of gene sets (each segment -> its single closest gene, then the
unique set of those genes within the segment's BMAP quartile). A gene's quartile is inherited from
the segment, so a gene can appear in more than one quartile. This is the canonical, extended
GO / KEGG / Reactome / Hallmark categories and for membrane glycoproteins.

Categories tested: GO_{BP,CC,MF}_2023, KEGG_2021_Human, Reactome_2022, MSigDB_Hallmark_2020
(from Enrichr GMTs), plus UniProt KW-0325 glycoproteins and membrane glycoproteins (KW-0325 &
KW-0472). Within each quartile only categories with ktest >= MINK are tested; BH-FDR across
categories within a quartile, separately for the hypergeometric and permutation p-values.
MHC excluded from both test and background.

Writes: results/uncoal_enrichment_both_Q{1..4}.tsv  and  results/uncoal_enrichment_glyco_summary.tsv
"""
import os
import glob
import numpy as np
import pandas as pd
from pathlib import Path

HERE = Path(__file__).resolve().parent
from scipy.stats import hypergeom

AP = str(HERE / "resources")
# gene-set libraries (.gmt) live with the gene-set enrichment
GA = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                  "../../6_TSP_Analyses/enrichment/geneset")
GLYCO = "../../6_TSP_Analyses/enrichment/glycoprotein/human_glycoproteins_KW0325_with_coords.tsv"
MEMGLYCO = "../../6_TSP_Analyses/enrichment/glycoprotein/human_membrane_glycoproteins_KW0325_KW0472_with_coords.tsv"
MHC = ("chr6", 28510120, 33480577)
DIST = None       # PRIMARY: no segment->gene distance filter 
NPERM = 1000
MINK = 3          # a category needs >= MINK observed test genes to be tested

def load_gmt(p):
    d = {}
    for line in open(p):
        f = line.rstrip("\n").split("\t")
        if len(f) > 2:
            d[f[0]] = {g.upper() for g in f[2:] if g}
    return d


def load_geneset_tsv(path):
    gl = pd.read_csv(path, sep="\t")
    col = [c for c in gl.columns
           if c.lower() in ("hgnc_symbol", "gene_name_primary", "symbol", "gene names (primary)")][0]
    return {str(x).upper() for x in gl[col].dropna()}


# gene-set libraries: Enrichr GMTs + the two UniProt keyword sets as a single "UniProt" library
libs = {g.split("/")[-1][:-4]: load_gmt(g) for g in sorted(glob.glob(f"{GA}/genesets/*.gmt"))}
libs["UniProt"] = {
    "KW-0325 glycoprotein": load_geneset_tsv(GLYCO),
    "membrane glycoprotein (KW-0325 & KW-0472)": load_geneset_tsv(MEMGLYCO),
}

# BMAP quartile edges from the length-weighted BMAP distribution of ALL decoded segments
bc = pd.read_csv(f"{AP}/all_segments/intersect/all_segments_intersect.whole_genome.bmap.length_counts.txt",
                 header=None, sep=" ", names=["bmap", "count"]).sort_values("bmap")
bc = bc[bc.bmap != "NA"].copy()
bc["bmap"] = bc.bmap.astype(float)
bc["cf"] = bc["count"].cumsum() / bc["count"].sum()
edges = [bc[bc.cf >= q].bmap.iloc[0] for q in [0.25, 0.5, 0.75]]


def quart(b):
    return 1 if b <= edges[0] else 2 if b <= edges[1] else 3 if b <= edges[2] else 4


def load_seg(path, cols):
    df = pd.read_csv(path, header=None, sep="\t")
    df.columns = cols
    df["distance"] = pd.to_numeric(df.distance, errors="coerce")
    if DIST is not None:
        df = df[df.distance <= DIST]
    df["bmap"] = pd.to_numeric(df.bmap, errors="coerce")
    df = df.dropna(subset=["bmap"])
    df = df[~((df.chrom == MHC[0]) & (df.start < MHC[2]) & (df.end > MHC[1]))]   # MHC excluded
    df["q"] = df.bmap.apply(quart)
    df["gene"] = df.gene.astype(str).str.upper()
    return df[["gene", "q"]]


test = load_seg(f"{AP}/expected_uncoal/intersect/expected_uncoal_intersect.whole_genome.exact.merged.closest_genes.bed",
                ["chrom", "start", "end", "p0", "phylop", "bmap", "cg", "gs", "ge", "gene", "distance"])
bg = load_seg(f"{AP}/all_segments/intersect/all_segments_intersect.whole_genome.closest_genes.dedup.bed",
              ["chrom", "start", "end", "phylop", "bmap", "ov", "cg", "gs", "ge", "gene", "distance"])


def bh(p):
    p = np.asarray(p, float)
    m = len(p)
    if m == 0:
        return p
    o = np.argsort(p)
    r = np.minimum.accumulate((p[o] * m / (np.arange(m) + 1))[::-1])[::-1]
    out = np.empty(m)
    out[o] = np.clip(r, 0, 1)
    return out


glyco_summary = []
for q in [1, 2, 3, 4]:
    tg = set(test[test.q == q].gene)
    nseg = len(test[test.q == q])            # # uncoalesced segments in quartile (permutation sample size)
    ntest = len(tg)                          # # unique test genes
    bgq = bg[bg.q == q]
    bg_arr = bgq.gene.values                 # one slot per background segment (its nearest gene)
    bg_genes = set(bgq.gene)
    nbg = len(bg_genes)                       # # unique background genes

    rng = np.random.default_rng(0)
    null_sets = [set(rng.choice(bg_arr, size=nseg, replace=False)) for _ in range(NPERM)]

    rows = []
    for lib, cats in libs.items():
        for cat, genes in cats.items():
            ktest = len(genes & tg)
            if ktest < MINK:
                continue
            kbg = len(genes & bg_genes)
            # (1) naive hypergeometric
            fold_hyperg = (ktest / ntest) / (kbg / nbg) if kbg else np.inf
            p_hyperg = float(hypergeom.sf(ktest - 1, nbg, kbg, ntest)) if kbg else 1.0
            # (2) segment permutation
            knull = np.array([len(genes & s) for s in null_sets])
            fold_perm = ktest / max(knull.mean(), 1e-9)
            perm_p = (1 + int((knull >= ktest).sum())) / (1 + NPERM)
            rows.append(dict(library=lib, category=cat, ktest=ktest, ntest=ntest, kbg=kbg, nbg=nbg,
                             fold_hyperg=round(fold_hyperg, 3), p_hyperg=p_hyperg,
                             null_mean=round(float(knull.mean()), 2),
                             fold_perm=round(fold_perm, 3), perm_p=round(perm_p, 4),
                             genes=",".join(sorted(genes & tg))))
    r = pd.DataFrame(rows).sort_values("perm_p").reset_index(drop=True)
    if len(r):
        r["p_hyperg_BH"] = bh(r["p_hyperg"])
        r["perm_p_BH"] = bh(r["perm_p"])
    r.to_csv(f"{AP}/results/uncoal_enrichment_both_Q{q}.tsv", sep="\t", index=False)

    for cat in ["KW-0325 glycoprotein", "membrane glycoprotein (KW-0325 & KW-0472)"]:
        row = r[r.category == cat]
        if len(row):
            row = row.iloc[0]
            glyco_summary.append(dict(quartile=f"Q{q}", category=cat, n_uncoal_segments=nseg,
                                      ktest=row.ktest, ntest=ntest, kbg=row.kbg, nbg=nbg,
                                      fold_hyperg=row.fold_hyperg, p_hyperg=round(row.p_hyperg, 4),
                                      fold_perm=row.fold_perm, perm_p=row.perm_p))
    print(f"Q{q}: {nseg} uncoal segs -> {ntest} test genes | {nbg} bg genes | "
          f"categories tested {len(r)}, perm BH<0.05 {int((r.perm_p_BH < 0.05).sum()) if len(r) else 0}",
          flush=True)

pd.DataFrame(glyco_summary).to_csv(f"{AP}/results/uncoal_enrichment_glyco_summary.tsv",
                                   sep="\t", index=False)
print(f"\nWrote results/uncoal_enrichment_both_Q1-4.tsv and uncoal_enrichment_glyco_summary.tsv")
