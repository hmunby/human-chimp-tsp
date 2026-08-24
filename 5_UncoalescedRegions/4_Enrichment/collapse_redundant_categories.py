#!/usr/bin/env python3
"""
Gene-set libraries are highly redundant, so a count of significant categories mostly measures how
many labels the library happens to attach to whichever genes you drew. This collapses categories
by the genes actually driving them, three ways:

  exact      : categories whose driver-gene sets are identical (misses heavy-but-imperfect
               overlap) 
               cardiac/membrane categories with slightly different extra genes)
  jaccard    : clustering at Jaccard >= T on driver-gene sets. Transitive, so a
               chain of partially-overlapping categories merges into one cluster.
  greedy     : minimum categories needed to cover every driver gene at least once (set cover).
               This is the number of independent findings.

Driver genes = category ∩ the observed gene set, restricted to genes retained in >=50% of the
tie-break resolutions.

Usage: python3 collapse_redundant_categories.py [D_kb ...]     (default 0 10)
"""

from pathlib import Path
HERE = Path(__file__).resolve().parent
# Repo-relative paths.
AP = str(HERE / "resources")
GA = str(HERE / "../../6_TSP_Analyses/enrichment/geneset")
GL = str(HERE / "../../6_TSP_Analyses/enrichment/glycoprotein")
import glob, sys
import numpy as np, pandas as pd


MHC=("chr6",28510120,33480577); NTIE=200
DS=[int(x) for x in sys.argv[1:]] or [0,10]

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
qv=lambda b: int(np.searchsorted(E,b,side="left")+1)
C=["chrom","start","end","p0","phylop","bmap","cg","gs","ge","gene","distance"]
d=pd.read_csv(f"{AP}/expected_uncoal/intersect/expected_uncoal_intersect.whole_genome.exact.merged.closest_genes.bed",
              header=None,sep="\t",names=C)
d["bmap"]=pd.to_numeric(d.bmap,errors="coerce"); d=d.dropna(subset=["bmap"])
d=d[~((d.chrom==MHC[0])&(d.start<MHC[2])&(d.end>MHC[1]))]
d["gene"]=d.gene.astype(str).str.upper()
reg=d.groupby(["chrom","start","end"],as_index=False).agg(bmap=("bmap","first")).sort_values(["chrom","start"]).reset_index(drop=True)
cand=d.groupby(["chrom","start","end"]).gene.apply(lambda s: sorted(set(s))).to_dict()
same=reg.chrom.values[1:]==reg.chrom.values[:-1]
gap=np.where(same, reg.start.values[1:]-reg.end.values[:-1], np.inf)

def freq(D):
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
    return out

def sim(a, b, how):
    """jaccard = |A&B|/|AuB|; overlap = |A&B|/min(|A|,|B|) (GOMCL, robust to size disparity);
    combined = 0.5*jaccard + 0.5*overlap (Enrichment Map default, threshold 0.375)."""
    i = len(a & b)
    if not i: return 0.0
    j = i / len(a | b)
    o = i / min(len(a), len(b))
    return {"jaccard": j, "overlap": o, "combined": 0.5*j + 0.5*o}[how]


def jaccard_clusters(sets, T, how="jaccard"):
    n=len(sets); par=list(range(n))
    def find(a):
        while par[a]!=a: par[a]=par[par[a]]; a=par[a]
        return a
    for i in range(n):
        for j in range(i+1,n):
            if sim(sets[i], sets[j], how) >= T:
                a,b=find(i),find(j)
                if a!=b: par[a]=b
    cl={}
    for i in range(n): cl.setdefault(find(i),[]).append(i)
    return list(cl.values())

def greedy_cover(sets):
    need=set().union(*sets) if sets else set()
    chosen=[]
    while need:
        best=max(range(len(sets)), key=lambda i: len(sets[i]&need))
        if not (sets[best]&need): break
        chosen.append(best); need-=sets[best]
    return chosen

SW=pd.concat([pd.read_csv(f,sep="\t") for f in
              [f"{AP}/results/uncoal_enrichment_locus_sweep.D0_D10k.tsv"]])
print(f"{'D':>4} {'Q':3s} {'cats':>5s} {'genes':>6s} | {'exact':>6s} {'J>=0.5':>7s} {'J>=0.3':>7s} {'Ovl>=.5':>7s} {'EM>=.375':>8s} {'greedy':>7s}")
detail={}
for D in DS:
    GF=freq(D*1000)
    for q in [1,2,3,4]:
        s=SW[(SW.D==D*1000)&(SW.quartile==q)&(SW.MIN_K==3)]
        s=s[(s.frac_pass>=.95)&(s.frac_tested>=.95)].sort_values("p_median")
        if not len(s): 
            print(f"{D:4d} Q{q:<2d} {0:5d} {0:6d} |      -       -       -       -        -       -"); continue
        f=GF[q]
        names=list(s.category); sets=[{g for g in CATS.get(c,set()) if f.get(g,0)>=NTIE*.5} for c in names]
        allg=set().union(*sets)
        exact=len({frozenset(x) for x in sets})
        c5=jaccard_clusters(sets,0.5); c3=jaccard_clusters(sets,0.3)
        ov=jaccard_clusters(sets,0.5,"overlap"); em=jaccard_clusters(sets,0.375,"combined")
        gc=greedy_cover(sets)
        print(f"{D:4d} Q{q:<2d} {len(s):5d} {len(allg):6d} | {exact:6d} {len(c5):7d} {len(c3):7d} {len(ov):7d} {len(em):7d} {len(gc):7d}")
        detail[(D,q)]=(names,sets,em,gc,s)
print("\n\n=== D=10 kb: Enrichment-Map clusters (0.5*Jaccard+0.5*Overlap >= 0.375) ===")
for q in [1,2,3,4]:
    if (10,q) not in detail: continue
    names,sets,c3,gc,s=detail[(10,q)]
    print(f"\n--- Q{q}: {len(names)} categories -> {len(c3)} clusters ---")
    for cl in sorted(c3,key=len,reverse=True):
        rep=cl[0]
        u=set().union(*[sets[i] for i in cl])
        print(f"  [{len(cl):2d} cat] {names[rep][:52]:52s} genes({len(u)}): {', '.join(sorted(u)[:10])}{' …' if len(u)>10 else ''}")
        for i in cl[1:6]: print(f"           + {names[i][:60]}")
        if len(cl)>6: print(f"           + {len(cl)-6} more")

# ── write the collapsed table (EnrichmentMap combined coefficient, k=0.5, threshold 0.375) ────
out=[]
for (D,q),(names,sets,cl,gc,srows) in sorted(detail.items()):
    srows=srows.reset_index(drop=True)
    for ci,members in enumerate(sorted(cl, key=len, reverse=True), 1):
        sub=srows.iloc[members]
        rep=sub.p_median.idxmin() if hasattr(sub.p_median,'idxmin') else members[0]
        rep_row=srows.loc[rep] if rep in srows.index else srows.iloc[members[0]]
        u=set().union(*[sets[i] for i in members])
        out.append(dict(
            merge_kb=D, quartile=f"Q{q}", cluster=ci, n_categories=len(members),
            representative=rep_row.category,
            best_BH=round(sub.BH_median.min(),4),
            fold_min=round(sub.fold.min(),2), fold_max=round(sub.fold.max(),2),
            k_min=int(sub.k_median.min()), k_max=int(sub.k_median.max()),
            n_genes=len(u), genes="; ".join(sorted(u)),
            all_categories=" | ".join(srows.iloc[i].category for i in members)))
T=pd.DataFrame(out)
fn=f"{AP}/results/uncoal_categories_collapsed_EM0.375.tsv"
T.to_csv(fn,sep="\t",index=False)
print(f"\n\nwrote {fn}  ({len(T)} clusters)")
for D in sorted(T.merge_kb.unique()):
    print(f"\n{'='*100}\nMERGE {D} kb, no distance cap — EnrichmentMap combined coefficient >= 0.375\n{'='*100}")
    for q in ["Q1","Q2","Q3","Q4"]:
        s2=T[(T.merge_kb==D)&(T.quartile==q)]
        if not len(s2): continue
        print(f"\n{q}: {s2.n_categories.sum()} categories -> {len(s2)} clusters")
        for _,x in s2.iterrows():
            print(f"  [{x.cluster}] {x.n_categories} cat | fold {x.fold_min}-{x.fold_max} | BH {x.best_BH} | {x.n_genes} genes")
            print(f"      {x.representative}")
            for c in x.all_categories.split(" | ")[1:]:
                print(f"      + {c}")
            print(f"      genes: {x.genes}")
