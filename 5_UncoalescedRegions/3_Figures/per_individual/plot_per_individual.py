#!/usr/bin/env python3
"""Per-individual uncoalesced-at-5.5-Mya by BMAP decile. 9 chimps coloured by subspecies
(EAS blue / WES purple / CEN red), 7 humans green.
Two modes:
  pct   -> y = % of genome uncoalesced (numerator/denominator*100)
  count -> y = expected number of 1 kb segments uncoalesced (numerator/1000)
Usage: plot_per_individual.py [DATADIR=data] [OUTSUFFIX=''] [MODE=pct]
Per-sample table columns: decile, denominator_bp, numerator(=sum posterior*overlap_bp), fraction."""
import sys, numpy as np, matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.cm as cm

BASE="../../../5_UncoalescedRegions/3_Figures/per_individual"
DATADIR   = sys.argv[1] if len(sys.argv) > 1 else "data"
OUTSUFFIX = sys.argv[2] if len(sys.argv) > 2 else ""
MODE      = sys.argv[3] if len(sys.argv) > 3 else "pct"     # 'pct' or 'count'

chimp = {"EAS_2003":"EAS","EAS_SAMEA4374790":"EAS","EAS_SAMEA4374797":"EAS",
         "WES_1059":"WES","WES_SAMEA5204228":"WES","WES_SAMEA2421542":"WES",
         "CEN_SAMEA4374772":"CEN","CEN_SAMEA4374778":"CEN","CEN_SAMEA4374785":"CEN"}
humans = ["ESN_HG03515","YRI_NA18488","MSL_HG03212","GWD_HG02568","ACB_HG01882","ASW_NA19625","LWK_NA19017"]

def shades(cmap,n): return [cm.get_cmap(cmap)(x) for x in np.linspace(0.55,0.9,n)]
col={}
for grp,cmap in [("EAS","Blues"),("WES","Purples"),("CEN","Reds")]:
    members=[s for s in chimp if chimp[s]==grp]
    for s,c in zip(members,shades(cmap,len(members))): col[s]=c
for s,c in zip(humans,shades("Greens",len(humans))): col[s]=c

def load(s):
    a=np.loadtxt(f"{BASE}/{DATADIR}/{s}.decile_frac.tsv")
    y = a[:,3]*100 if MODE=="pct" else a[:,2]/1000.0   # % of genome  OR  expected # 1kb segments
    return a[:,0], y

fig,ax=plt.subplots(figsize=(11,5.2),dpi=400)
# connect each individual's decile points with a line so overlapping markers are separable
for s in humans:
    d,y=load(s); ax.plot(d, y, '-^', color=col[s], ms=7, mec='white', mew=0.4, lw=1.2, alpha=0.85)
for s in chimp:
    d,y=load(s); ax.plot(d, y, '-o', color=col[s], ms=8, mec='white', mew=0.4, lw=1.2, alpha=0.85)

ax.set_xlabel("Background selection decile", fontweight="bold")
ax.set_ylabel("% of genome uncoalesced at 5.5 Mya" if MODE=="pct"
              else "Expected number of segments uncoalesced at 5.5 Mya", fontweight="bold")
ax.set_xticks(range(1,11))
_xtl=[str(b) for b in range(1,11)]
_xtl[0]="1\n"+r"$\mathbf{(strongest)}$"; _xtl[-1]="10\n"+r"$\mathbf{(weakest)}$"
ax.set_xticklabels(_xtl)
ax.set_xlim(0.5,10.5); ax.set_ylim(bottom=0)
leg=[Line2D([0],[0],marker='^',color=cm.get_cmap("Greens")(0.7),markerfacecolor=cm.get_cmap("Greens")(0.7),ms=10,lw=1.2,label='Human'),
     Line2D([0],[0],marker='o',color=cm.get_cmap("Blues")(0.7),markerfacecolor=cm.get_cmap("Blues")(0.7),ms=10,lw=1.2,label='Chimp (EAS)'),
     Line2D([0],[0],marker='o',color=cm.get_cmap("Purples")(0.7),markerfacecolor=cm.get_cmap("Purples")(0.7),ms=10,lw=1.2,label='Chimp (WES)'),
     Line2D([0],[0],marker='o',color=cm.get_cmap("Reds")(0.7),markerfacecolor=cm.get_cmap("Reds")(0.7),ms=10,lw=1.2,label='Chimp (CEN)')]
ax.legend(handles=leg, loc='upper left', frameon=False)
plt.tight_layout()
out=f"{BASE}/figures/SuppFig_per_individual_uncoal_bmap_decile{OUTSUFFIX}.png"
plt.savefig(out,dpi=400,bbox_inches="tight"); print("saved",out)
