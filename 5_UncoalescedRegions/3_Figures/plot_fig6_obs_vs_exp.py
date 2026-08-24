#!/usr/bin/env python3
"""Figure 6: observed vs expected reciprocal uncoalescence by background-selection quartile,
and the enrichment (observed/expected). MHC excluded.
"""
import os, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt

AP="resources"
os.chdir(AP)
cc="0.331998"; n_bins=4
LAB=["Q1\n(strongest)","Q2","Q3","Q4\n(weakest)"]
XLABEL="Background selection quartile"

def read_freqs(p):
    d={}
    for line in open(p):
        c,b=line.strip().split(); d[b]=int(c)
    return d
human_all=read_freqs(f"all_segments/human/phylop/all_segments_human_unconstrained_cutoff_{cc}_phylop100way.bmap.counts.txt")

def calc_bins(freqs,n):
    vals=[]
    for b,c in freqs.items():
        if b!="NA": vals+=[float(b)]*c
    return np.percentile(vals, np.linspace(0,100,n+1))
edges=calc_bins(human_all,n_bins)

def bin_len_from_freqs(edges,freqs,seglen):
    out={}
    for i in range(len(edges)-1):
        lo,hi=edges[i],edges[i+1]
        out[f"Bin {i+1}"]=sum(c for b,c in freqs.items() if b!="NA" and lo<=float(b)<hi)*seglen
    return out
human_bin_lengths=bin_len_from_freqs(edges,human_all,1000)

def bin_len_from_lc(p,edges):
    out={f"Bin {i+1}":0 for i in range(len(edges)-1)}
    for line in open(p):
        b,l=line.strip().split()
        if b=="NA": continue
        bv=float(b)
        for i in range(len(edges)-1):
            if edges[i]<=bv<edges[i+1]: out[f"Bin {i+1}"]+=int(l); break
    return out
chimp_bin_lengths=bin_len_from_lc(f"all_segments/chimp/phylop/all_segments_chimp_unconstrained_length_cutoff_{cc}_phylop100way.bmap.length_counts.txt",edges)

labels=[f"Bin {i+1}" for i in range(n_bins)]
def cut(x): return pd.cut(x,bins=edges,labels=labels,include_lowest=True)

chimp_df=pd.read_csv(f"expected_uncoal/chimp/expected_uncoal_chimp_5755625.genome_wide.phylop100way.unconstrained_cutoff_{cc}.bed",sep="\t",header=None,names=["chrom","start","end","p0","phylop","bmap"])
chimp_df["len"]=chimp_df.end-chimp_df.start
human_df=pd.read_csv(f"expected_uncoal/human/phylop/expected_uncoal_human_5589889.genome_wide.phylop100way.unconstrained_cutoff_{cc}.bed",sep="\t",header=None,names=["chrom","start","end","p0","phylop","bmap"])
human_df["len"]=human_df.end-human_df.start
inter=pd.read_csv(f"expected_uncoal/intersect/expected_uncoal_intersect.unconstrained_cutoff_{cc}.bed",sep="\t",header=None,
    names=["chrom_h","start_h","end_h","p0_h","ph_h","bmap_h","chrom_c","start_c","end_c","id_c","ph_c","bmap_c","overlap_len"])
inter=inter[inter.overlap_len>0].reset_index(drop=True)

def rm_mhc(df,chrom,start,end):
    return df[~((df[chrom]=="chr6")&(df[start]<34000000)&(df[end]>28000000))].copy()
chimp_m=rm_mhc(chimp_df,"chrom","start","end"); human_m=rm_mhc(human_df,"chrom","start","end")
inter_m=rm_mhc(inter,"chrom_h","start_h","end_h")
chimp_m["bb"]=cut(chimp_m["bmap"]); human_m["bb"]=cut(human_m["bmap"]); inter_m["bb"]=cut(inter_m["bmap_h"])

def fracs(chimp_d,human_d,inter_d):
    cl=chimp_d.groupby("bb")["len"].sum().to_dict(); hl=human_d.groupby("bb")["len"].sum().to_dict()
    ov=inter_d.groupby("bb")["overlap_len"].sum().to_dict()
    cf={b: cl.get(b,0)/chimp_bin_lengths[b] if chimp_bin_lengths[b]>0 else 0 for b in labels}
    hf={b: hl.get(b,0)/human_bin_lengths[b] if human_bin_lengths[b]>0 else 0 for b in labels}
    of={b: ov.get(b,0)/chimp_bin_lengths[b] if chimp_bin_lengths[b]>0 else 0 for b in labels}
    ef={b: hf[b]*cf[b] for b in labels}
    en={b: (of[b]/ef[b]) if ef[b]>0 else np.nan for b in labels}
    return of,ef,en
obs,exp,enr=fracs(chimp_m,human_m,inter_m)
print("no-MHC quartile enrichment (obs/exp):",{b:round(enr[b],2) for b in labels})

def ci(d): return np.nanpercentile(d,2.5),np.nanpercentile(d,97.5)
nb=200

# chromosome-block resampling: pick whole chromosomes with replacement and apply the SAME
# set to chimp/human/intersect, preserving numerator-denominator coupling and within-chromosome
# autocorrelation. (Segment resampling instead treats each segment as independent.)
_grp={"chimp":{c:g for c,g in chimp_m.groupby("chrom")},
      "human":{c:g for c,g in human_m.groupby("chrom")},
      "inter":{c:g for c,g in inter_m.groupby("chrom_h")}}
_chroms=sorted(set(_grp["human"])|set(_grp["chimp"])|set(_grp["inter"]))

def _resample(method,rng):
    if method=="segment":
        return (chimp_m.sample(frac=1,replace=True,random_state=rng),
                human_m.sample(frac=1,replace=True,random_state=rng),
                inter_m.sample(frac=1,replace=True,random_state=rng))
    picks=rng.choice(_chroms,size=len(_chroms),replace=True)
    def cat(g):
        parts=[g[c] for c in picks if c in g]
        return pd.concat(parts) if parts else next(iter(g.values())).iloc[:0]
    return cat(_grp["chimp"]),cat(_grp["human"]),cat(_grp["inter"])

def _bootstrap(method):
    ob={b:[] for b in labels}; ex={b:[] for b in labels}; en={b:[] for b in labels}
    rng=np.random.RandomState(0)
    for _ in range(nb):
        o,e,en2=fracs(*_resample(method,rng))
        for b in labels:
            ob[b].append(o[b]); ex[b].append(e[b]); en[b].append(en2[b])
    return ob,ex,en

def build(ob_bs,ex_bs,en_bs,out):
    x=np.arange(1,n_bins+1)
    fig,(axA,axB)=plt.subplots(1,2,figsize=(11,4.6),dpi=200)
    o=[obs[b]*100 for b in labels]; e=[exp[b]*100 for b in labels]
    oerr=np.array([[(obs[b]-ci(ob_bs[b])[0])*100 for b in labels],[(ci(ob_bs[b])[1]-obs[b])*100 for b in labels]])
    eerr=np.array([[(exp[b]-ci(ex_bs[b])[0])*100 for b in labels],[(ci(ex_bs[b])[1]-exp[b])*100 for b in labels]])
    axA.errorbar(x,o,yerr=oerr,fmt="-o",color="#1F77B4",capsize=3,label="Observed")
    axA.errorbar(x,e,yerr=eerr,fmt="-o",color="#7F7F7F",capsize=3,label="Expected")
    axA.set_xlabel(XLABEL,fontweight="bold")
    axA.set_ylabel("Percentage of genome uncoalesced in both species",fontweight="bold")
    axA.set_xticks(x); axA.set_xticklabels(LAB); axA.legend(frameon=False)
    axA.set_title("A",loc="left",fontweight="bold")
    axA.spines[["top","right"]].set_visible(False)
    en_v=[enr[b] for b in labels]
    enerr=np.array([[enr[b]-ci(en_bs[b])[0] for b in labels],[ci(en_bs[b])[1]-enr[b] for b in labels]])
    axB.errorbar(x,en_v,yerr=enerr,fmt="-o",color="#9467BD",capsize=3)
    axB.axhline(1,ls="--",color="grey",lw=1)
    axB.set_xlabel(XLABEL,fontweight="bold")
    axB.set_ylabel("Enrichment (observed / expected)",fontweight="bold")
    axB.set_xticks(x); axB.set_xticklabels(LAB); axB.set_ylim(0,max(en_v)*1.3)
    axB.set_title("B",loc="left",fontweight="bold")
    axB.spines[["top","right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out,dpi=200,bbox_inches="tight")
    print("saved",out)

# segment bootstrap -> main figure; chromosome-block -> separate ".chromblock" figure
for _method,_suffix in [("segment",""),("chromosome",".chromblock")]:
    _ob,_ex,_en=_bootstrap(_method)
    print(f"{_method} enrichment 95% CI:",{b:(round(ci(_en[b])[0],2),round(ci(_en[b])[1],2)) for b in labels})
    build(_ob,_ex,_en,f"{AP}/reciprocal_uncoal_obs_vs_exp.quartiles{_suffix}.png")
