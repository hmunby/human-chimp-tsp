#!/usr/bin/env python3
"""
Gene-level ("naive") hypergeometric for the two glycoprotein sets, across every merge distance and
distance cap.

Test set  : unique genes nearest to the uncoalesced LOCI in a BMAP quartile, one gene per locus,
            averaged over N_TIEBREAK random resolutions of bedtools-closest ties.
Background: unique genes nearest to ANY callable segment in the same quartile (MHC excluded, and
            the same distance cap applied).
fold = (ktest/ntest) / (kbg/nbg);  p = one-sided hypergeometric P(X >= ktest).
"""

from pathlib import Path
HERE = Path(__file__).resolve().parent
# Repo-relative paths. Large intermediates (the uncoalesced-segment and all-segment BEDs, ~4.9M
# rows) are NOT bundled; place them under resources/.
AP = str(HERE / "resources")
GA = str(HERE / "../../6_TSP_Analyses/enrichment/geneset")
GL = str(HERE / "../../6_TSP_Analyses/enrichment/glycoprotein")
import numpy as np, pandas as pd
from scipy.stats import hypergeom


MHC=("chr6",28510120,33480577); NTIE=200 

def tsvset(p):
    g=pd.read_csv(p,sep="\t"); c=[x for x in g.columns if x.lower() in
        ("hgnc_symbol","gene_name_primary","symbol","gene names (primary)")][0]
    return {str(x).upper() for x in g[c].dropna()}
SETS={"KW-0325 glycoprotein":tsvset(f"{GL}/human_glycoproteins_KW0325_with_coords.tsv"),
      "membrane glycoprotein (KW-0325 & KW-0472)":tsvset(f"{GL}/human_membrane_glycoproteins_KW0325_KW0472_with_coords.tsv")}

bc=pd.read_csv(f"{AP}/all_segments/intersect/all_segments_intersect.whole_genome.bmap.length_counts.txt",
               header=None,sep=" ",names=["bmap","count"]).sort_values("bmap")
bc=bc[bc.bmap!="NA"].copy(); bc["bmap"]=bc.bmap.astype(float); bc["cf"]=bc["count"].cumsum()/bc["count"].sum()
E=np.array([bc[bc.cf>=q].bmap.iloc[0] for q in [.25,.5,.75]])
qv=lambda b: np.searchsorted(E,np.asarray(b,float),side="left")+1

C=["chrom","start","end","p0","phylop","bmap","cg","gs","ge","gene","distance"]
T=pd.read_csv(f"{AP}/expected_uncoal/intersect/expected_uncoal_intersect.whole_genome.exact.merged.closest_genes.bed",
              header=None,sep="\t",names=C)
T["bmap"]=pd.to_numeric(T.bmap,errors="coerce"); T=T.dropna(subset=["bmap"])
T=T[~((T.chrom==MHC[0])&(T.start<MHC[2])&(T.end>MHC[1]))]
T["gene"]=T.gene.astype(str).str.upper(); T["distance"]=pd.to_numeric(T.distance,errors="coerce")

BC=["chrom","start","end","phylop","bmap","ov","cg","gs","ge","gene","distance"]
B=pd.read_csv(f"{AP}/all_segments/intersect/all_segments_intersect.whole_genome.closest_genes.dedup.bed",
              header=None,sep="\t",names=BC,usecols=["chrom","start","end","bmap","gene","distance"])
B["bmap"]=pd.to_numeric(B.bmap,errors="coerce"); B=B.dropna(subset=["bmap"])
B=B[~((B.chrom==MHC[0])&(B.start<MHC[2])&(B.end>MHC[1]))]
B["gene"]=B.gene.astype(str).str.upper(); B["distance"]=pd.to_numeric(B.distance,errors="coerce").fillna(1<<30)
B["q"]=qv(B.bmap.values)

rows=[]
for cap in [np.inf,100_000,50_000]:
    t = T if not np.isfinite(cap) else T[T.distance<=cap]
    b = B if not np.isfinite(cap) else B[B.distance<=cap]
    bgq={q:set(b[b.q==q].gene) for q in [1,2,3,4]}
    reg=t.groupby(["chrom","start","end"],as_index=False).agg(bmap=("bmap","first")).sort_values(["chrom","start"]).reset_index(drop=True)
    cand=t.groupby(["chrom","start","end"]).gene.apply(lambda s: sorted(set(s))).to_dict()
    same=reg.chrom.values[1:]==reg.chrom.values[:-1]
    gap=np.where(same, reg.start.values[1:]-reg.end.values[:-1], np.inf)
    for D in [0,5000,10000,25000]:
        newl=np.ones(len(reg),dtype=bool); newl[1:]=(~same)|(gap>D)
        reg["lid"]=np.cumsum(newl)-1
        loci=[]
        for _,g in reg.groupby("lid"):
            w=(g.end-g.start).values.astype(float)
            cs=sorted({x for _,r in g.iterrows() for x in cand[(r.chrom,r.start,r.end)]})
            loci.append((int(qv(np.average(g.bmap.values,weights=w))),cs))
        for q in [1,2,3,4]:
            cs=[c for qq,c in loci if qq==q]
            if not cs: continue
            rng=np.random.default_rng(500+q)
            sets=[{c[rng.integers(0,len(c))] for c in cs} for _ in range(NTIE)]
            nbg=len(bgq[q])
            for name,gs in SETS.items():
                kb=len(gs&bgq[q])
                kt=np.array([len(gs&s) for s in sets]); nt=np.array([len(s) for s in sets])
                fold=np.median((kt/nt)/(kb/nbg)); p=float(np.median([hypergeom.sf(a-1,nbg,kb,n) for a,n in zip(kt,nt)]))
                rows.append(dict(merge_kb=D//1000, dist_cap_kb=("none" if not np.isfinite(cap) else int(cap)//1000),
                                 quartile=f"Q{q}", category=name, n_loci=len(cs),
                                 ktest=int(np.median(kt)), ntest=int(np.median(nt)), kbg=kb, nbg=nbg,
                                 obs_pct=round(100*np.median(kt/nt),1), bg_pct=round(100*kb/nbg,1),
                                 fold_hyperg=round(fold,3), p_hyperg=p))
r=pd.DataFrame(rows)
fn=f"{AP}/results/glyco_genelevel_hypergeometric.tsv"
r.to_csv(fn,sep="\t",index=False)
print(f"wrote {fn}\n")
for cap in r.dist_cap_kb.unique():
    print(f"===== distance cap: {cap} =====")
    s=r[r.dist_cap_kb==cap]
    print(s.pivot_table(index=["quartile","category"],columns="merge_kb",
                        values=["fold_hyperg"],aggfunc="first").to_string())
    print()
