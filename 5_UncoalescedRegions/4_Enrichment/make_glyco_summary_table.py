#!/usr/bin/env python3
"""
Consolidated glycoprotein / membrane-glycoprotein enrichment table.

Joins the two tests for the same conditions:
  gene-level  = hypergeometric, uncoalesced loci vs all genes nearest to any callable segment in
                the same BMAP quartile (from glyco_genelevel_hypergeometric.tsv)
  permutation = span-matched locus null, FRACTION statistic, ties averaged over 200 resolutions
                (from the uncoal_enrichment_locus_sweep.* files)

Two BH columns:
  BH_scan   = Benjamini-Hochberg across every category tested in that quartile (~300-2500 tests).
              This is what the category scan reports. Conservative for these two sets, which were
              pre-specified rather than discovered by scanning.
  BH_glyco8 = Benjamini-Hochberg across just the 8 pre-specified tests reported here
              (2 gene sets x 4 quartiles), computed separately for each test and condition.
"""

from pathlib import Path
HERE = Path(__file__).resolve().parent
# Repo-relative paths.
GA = str(HERE / "../../6_TSP_Analyses/enrichment/geneset")
GL = str(HERE / "../../6_TSP_Analyses/enrichment/glycoprotein")
import glob, os, re
import numpy as np, pandas as pd


GLY=["KW-0325 glycoprotein","membrane glycoprotein (KW-0325 & KW-0472)"]

def bh(p):
    p=np.asarray(p,float); m=len(p)
    if m==0: return p
    o=np.argsort(p); r=np.minimum.accumulate((p[o]*m/(np.arange(m)+1))[::-1])[::-1]
    out=np.empty(m); out[o]=np.clip(r,0,1); return out

H=pd.read_csv(f"{AP}/results/glyco_genelevel_hypergeometric.tsv",sep="\t")
rows=[]
for f in sorted(glob.glob(f"{AP}/results/uncoal_enrichment_locus_sweep.*.tsv")):
    b=os.path.basename(f)
    m=re.search(r"cap(\d+)kb",b); cap="none" if not m else int(m.group(1))
    s=pd.read_csv(f,sep="\t"); s=s[(s.MIN_K==3)&(s.category.isin(GLY))]
    for _,x in s.iterrows():
        rows.append(dict(merge_kb=int(x.D)//1000, dist_cap_kb=cap, quartile=f"Q{int(x.quartile)}",
                         category=x.category, perm_fold=round(x.fold,3), perm_p=x.p_median,
                         perm_BH_scan=round(x.BH_median,4)))
P=pd.DataFrame(rows).drop_duplicates(subset=["merge_kb","dist_cap_kb","quartile","category"])
H["dist_cap_kb"]=H.dist_cap_kb.astype(str)
P["dist_cap_kb"]=P.dist_cap_kb.astype(str)
T=H.merge(P,on=["merge_kb","dist_cap_kb","quartile","category"],how="left")

# BH across the 8 pre-specified glyco tests, within each condition
for col,new in [("p_hyperg","hyperg_BH_glyco8"),("perm_p","perm_BH_glyco8")]:
    T[new]=np.nan
    for (mk,cap),g in T.groupby(["merge_kb","dist_cap_kb"]):
        v=g[col]
        if v.notna().sum()==0: continue
        idx=v.dropna().index
        T.loc[idx,new]=bh(T.loc[idx,col].values)
# no rounding: BH values here span 1e-8 to 1, rounding to 4dp turned the smallest into 0.0000

# labels make the nesting explicit: the membrane set (3,374 genes) is a strict SUBSET
# of all glycoproteins (4,531).
T["set"]=np.where(T.category.str.startswith("KW"),"glyco-all","glyco-membrane")
cols=["merge_kb","dist_cap_kb","quartile","set","category","n_loci","ktest","ntest","obs_pct","kbg","nbg","bg_pct",
      "fold_hyperg","p_hyperg","hyperg_BH_glyco8","perm_fold","perm_p","perm_BH_glyco8","perm_BH_scan"]
T=T[cols].sort_values(["dist_cap_kb","merge_kb","quartile","category"])
fn=f"{AP}/results/glyco_enrichment_SUMMARY.tsv"
T.to_csv(fn,sep="\t",index=False)
print(f"wrote {fn}  ({len(T)} rows)\n")
def show(mk,cap,title):
    s=T[(T.merge_kb==mk)&(T.dist_cap_kb==str(cap))]
    if not len(s): return
    print(f"\n{title}")
    print(f"{'Q':3s} {'set':15s} {'obs':>12s} {'bg':>7s} | {'GENE-LEVEL (hypergeometric)':^34s} | {'PERMUTATION (locus null, fraction)':^36s}")
    print(f"{'':3s} {'':9s} {'':12s} {'':7s} | {'fold':>6s} {'p':>10s} {'BH8':>7s} {'':7s} | {'fold':>6s} {'p':>10s} {'BH8':>7s} {'BHscan':>8s}")
    for _,x in s.iterrows():
        c=x.set
        pf="" if pd.isna(x.perm_fold) else f"{x.perm_fold:6.2f}"
        pp="" if pd.isna(x.perm_p) else f"{x.perm_p:10.2e}"
        p8="" if pd.isna(x.perm_BH_glyco8) else f"{x.perm_BH_glyco8:7.4f}"
        ps="" if pd.isna(x.perm_BH_scan) else f"{x.perm_BH_scan:8.4f}"
        print(f"{x.quartile:3s} {c:15s} {x.ktest:4d}/{x.ntest:<4d}={x.obs_pct:4.1f}% {x.bg_pct:6.1f}% | "
              f"{x.fold_hyperg:6.2f} {x.p_hyperg:10.2e} {x.hyperg_BH_glyco8:7.4f} {'':7s} | {pf} {pp} {p8} {ps}")
show(10,"none","PRIMARY: merge 10 kb, no distance cap")
show(0,"none","merge 0 kb (regions, not loci), no cap")
show(25,"none","merge 25 kb, no cap")
show(10,100,"merge 10 kb, 100 kb distance cap")
show(10,50,"merge 10 kb, 50 kb distance cap")
