#!/usr/bin/env python3
"""
One full-results TSV per merge distance: every significant category, its driver genes, and the
gene-overlap cluster it belongs to.

A category is significant when it passes BH<0.05 in >=95% of the 200 tie-break resolutions AND is
tested (reaches MIN_K) in >=95% of them. Both conditions are needed.

Clusters = connected components under the EnrichmentMap combined coefficient
(0.5*Jaccard + 0.5*overlap >= 0.375) on the driver genes (Merico et al. 2010, PLoS ONE 5:e13984).
Clusters are numbered per quartile per MIN_K, largest first; singletons get their own number.

Genes listed are those retained in >=50% of the tie-break resolutions; a percentage is appended to
any gene not retained in all of them.

Usage: python3 make_full_results_tables.py OUTDIR sweep1.tsv [sweep2.tsv ...]
"""

from pathlib import Path
HERE = Path(__file__).resolve().parent
# Repo-relative paths.
AP = str(HERE / "resources")
GA = str(HERE / "../../6_TSP_Analyses/enrichment/geneset")
GL = str(HERE / "../../6_TSP_Analyses/enrichment/glycoprotein")
import glob, os, re, sys
import numpy as np, pandas as pd


MHC=("chr6",28510120,33480577); NTIE=200; T=0.375
OUTDIR=sys.argv[1]; SWEEPS=sys.argv[2:]
os.makedirs(OUTDIR,exist_ok=True)

def gmt(p):
    o={}
    for l in open(p):
        f=l.rstrip("\n").split("\t")
        if len(f)>2: o[f[0]]={g.upper() for g in f[2:] if g}
    return o
def tsvset(p):
    g=pd.read_csv(p,sep="\t"); c=[x for x in g.columns if x.lower() in
        ("hgnc_symbol","gene_name_primary","symbol","gene names (primary)")][0]
    return {str(x).upper() for x in g[c].dropna()}
CATS={}
for f in sorted(glob.glob(f"{GA}/genesets/*.gmt")): CATS.update(gmt(f))
CATS["KW-0325 glycoprotein"]=tsvset(f"{GL}/human_glycoproteins_KW0325_with_coords.tsv")
CATS["membrane glycoprotein (KW-0325 & KW-0472)"]=tsvset(f"{GL}/human_membrane_glycoproteins_KW0325_KW0472_with_coords.tsv")

bc=pd.read_csv(f"{AP}/all_segments/intersect/all_segments_intersect.whole_genome.bmap.length_counts.txt",
               header=None,sep=" ",names=["bmap","count"]).sort_values("bmap")
bc=bc[bc.bmap!="NA"].copy(); bc["bmap"]=bc.bmap.astype(float); bc["cf"]=bc["count"].cumsum()/bc["count"].sum()
E=np.array([bc[bc.cf>=q].bmap.iloc[0] for q in [.25,.5,.75]])
qv=lambda b:int(np.searchsorted(E,b,side="left")+1)
C=["chrom","start","end","p0","phylop","bmap","cg","gs","ge","gene","distance"]
RAW=pd.read_csv(f"{AP}/expected_uncoal/intersect/expected_uncoal_intersect.whole_genome.exact.merged.closest_genes.bed",
                header=None,sep="\t",names=C)
RAW["bmap"]=pd.to_numeric(RAW.bmap,errors="coerce"); RAW=RAW.dropna(subset=["bmap"])
RAW=RAW[~((RAW.chrom==MHC[0])&(RAW.start<MHC[2])&(RAW.end>MHC[1]))]
RAW["gene"]=RAW.gene.astype(str).str.upper()
RAW["distance"]=pd.to_numeric(RAW.distance,errors="coerce")

_c={}
def gene_freq(D,cap=np.inf):
    if (D,cap) in _c: return _c[(D,cap)]
    d=RAW if not np.isfinite(cap) else RAW[RAW.distance<=cap]
    reg=d.groupby(["chrom","start","end"],as_index=False).agg(bmap=("bmap","first")).sort_values(["chrom","start"]).reset_index(drop=True)
    cand=d.groupby(["chrom","start","end"]).gene.apply(lambda s:sorted(set(s))).to_dict()
    same=reg.chrom.values[1:]==reg.chrom.values[:-1]
    gap=np.where(same,reg.start.values[1:]-reg.end.values[:-1],np.inf)
    newl=np.ones(len(reg),dtype=bool); newl[1:]=(~same)|(gap>D)
    reg["lid"]=np.cumsum(newl)-1
    loci=[]
    for _,g in reg.groupby("lid"):
        w=(g.end-g.start).values.astype(float)
        cs=sorted({x for _,r in g.iterrows() for x in cand[(r.chrom,r.start,r.end)]})
        loci.append((qv(np.average(g.bmap.values,weights=w)),cs))
    out={}
    for q in [1,2,3,4]:
        cs=[c for qq,c in loci if qq==q]; rng=np.random.default_rng(500+q); f={}
        for _ in range(NTIE):
            for g_ in {c[rng.integers(0,len(c))] for c in cs}: f[g_]=f.get(g_,0)+1
        out[q]=f
    _c[(D,cap)]=out; return out

def comb(a,b):
    i=len(a&b)
    if not i: return 0.0
    return 0.5*(i/len(a|b))+0.5*(i/min(len(a),len(b)))
def clusters(sets):
    n=len(sets); par=list(range(n))
    def f(a):
        while par[a]!=a: par[a]=par[par[a]]; a=par[a]
        return a
    for i in range(n):
        for j in range(i+1,n):
            if comb(sets[i],sets[j])>=T:
                a,b=f(i),f(j)
                if a!=b: par[a]=b
    g={}
    for i in range(n): g.setdefault(f(i),[]).append(i)
    return sorted(g.values(),key=lambda c:(-len(c),c[0]))

for path in SWEEPS:
    m=re.search(r"cap(\d+)kb",os.path.basename(path))
    cap=float(m.group(1))*1000 if m else np.inf
    sw=pd.read_csv(path,sep="\t")
    for D in sorted(sw.D.unique()):
        GF=gene_freq(int(D),cap); rows=[]
        for q in [1,2,3,4]:
            f=GF[q]
            for mk in sorted(sw.MIN_K.unique()):
                r=sw[(sw.D==D)&(sw.quartile==q)&(sw.MIN_K==mk)]
                r=r[(r.frac_pass>=.95)&(r.frac_tested>=.95)].sort_values("p_median").reset_index(drop=True)
                if not len(r): continue
                sets=[{g for g in CATS.get(c,set()) if f.get(g,0)>=NTIE*.5} for c in r.category]
                cl=clusters(sets)
                for ci,mem in enumerate(cl,1):
                    cu=set().union(*[sets[i] for i in mem])
                    for i in sorted(mem,key=lambda i:r.p_median[i]):
                        drv=sorted(sets[i],key=lambda g:(-f[g],g))
                        rows.append(dict(
                            merge_kb=int(D)//1000,
                            dist_cap_kb=("none" if not np.isfinite(cap) else int(cap)//1000),
                            quartile=f"Q{q}", MIN_K=int(mk),
                            cluster=ci, cluster_size=len(mem),
                            category=r.category[i], k=int(round(r.k_median[i])),
                            fold=round(r.fold[i],2), p=r.p_median[i],
                            BH=round(r.BH_median[i],4),
                            pct_resolutions_passing=int(round(100*r.frac_pass[i])),
                            n_genes=len(drv),
                            genes=", ".join(g if f[g]==NTIE else f"{g}({100*f[g]//NTIE}%)" for g in drv),
                            cluster_n_genes=len(cu),
                            cluster_genes="; ".join(sorted(cu))))
        t=pd.DataFrame(rows)
        tag=f"D{int(D)//1000}kb"+("" if not np.isfinite(cap) else f".cap{int(cap)//1000}kb")
        fn=f"{OUTDIR}/uncoal_full_results.{tag}.tsv"
        t.to_csv(fn,sep="\t",index=False)
        n3=t[t.MIN_K==3] if len(t) else t
        summ=", ".join(f"Q{q} {len(n3[n3.quartile==f'Q{q}'])} cats/{n3[n3.quartile==f'Q{q}'].cluster.nunique()} clusters"
                       for q in [1,2,3,4]) if len(t) else "empty"
        print(f"{fn}\n    MIN_K=3: {summ}", flush=True)
