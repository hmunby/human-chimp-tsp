#!/usr/bin/env python3
"""Parse Bitarello et al. 2018 supplementary tables S4 (candidate genes) and S5 (candidate windows)
into clean TSVs for the overlap comparison.

S5_Table_GBE.xlsx: 24 window sheets, one per {SIGNIF,OUTLIERS} x tf{0.3,0.4,0.5} x pop{LWK,YRI,GBR,TSI}.
  Each sheet is a headerless 3-col bed (chrom, start, end; hg19; 3 kb NCD windows).
S4_Table_GBE_revised.xlsx: candidate protein-coding genes overlapped by candidate windows, by
  continent x set: {Afr,Eur,Shared} x {outliers,significant}. Cols: Chr, Acronym, tf/p per population.

Writes: bitarello_windows.hg19.tsv  (chrom,start,end,set,tf,pop)
        bitarello_candidate_genes.tsv (chr,gene,set,continent,min_p)
"""
import re
import warnings
import pandas as pd

warnings.filterwarnings("ignore")
A = "../../../6_TSP_Analyses/literature/bitarello_2018"
POPS = ["LWK", "YRI", "GBR", "TSI"]

# --- S5 windows ------------------------------------------------------------------------------
xl5 = pd.ExcelFile(f"{A}/S5_Table_GBE.xlsx")
win = []
for sh in xl5.sheet_names:
    if not (sh.startswith("SIGNIF") or sh.startswith("OUTLIERS")):
        continue
    setname = "OUTLIERS" if sh.startswith("OUTLIERS") else "SIGNIF"
    tf = re.search(r"0\.+([345])", sh).group(1)          # handles 'tf0.5', 'tf0..3'
    pop = next(p for p in POPS if p in sh)
    d = pd.read_excel(xl5, sheet_name=sh, header=None).iloc[:, :3]
    d.columns = ["chrom", "start", "end"]
    d = d.dropna()
    d["set"], d["tf"], d["pop"] = setname, f"0.{tf}", pop
    win.append(d)
W = pd.concat(win, ignore_index=True)
W["chrom"] = W.chrom.astype(int); W["start"] = W.start.astype(int); W["end"] = W.end.astype(int)
W.to_csv(f"{A}/bitarello_windows.hg19.tsv", sep="\t", index=False)

# --- S4 candidate genes ----------------------------------------------------------------------
xl4 = pd.ExcelFile(f"{A}/S4_Table_GBE_revised.xlsx")
gene_sheets = {"Afr_outliers": ("Afr", "outliers"), "Eur_outliers": ("Eur", "outliers"),
               "Shared_outliers": ("Shared", "outliers"), "Afr_significant": ("Afr", "significant"),
               "Eur_significant": ("Eur", "significant"), "Shared_significant": ("Shared", "significant")}
genes = []
for sh, (cont, setname) in gene_sheets.items():
    d = pd.read_excel(xl4, sheet_name=sh)
    d = d[[c for c in d.columns if not str(c).startswith("Unnamed")]]
    pcols = [c for c in d.columns if str(c).startswith("p.")]
    d["min_p"] = d[pcols].min(axis=1)
    out = d[["Chr", "Acronym", "min_p"]].dropna(subset=["Acronym"]).copy()
    out["set"], out["continent"] = setname, cont
    genes.append(out)
G = pd.concat(genes, ignore_index=True).rename(columns={"Chr": "chr", "Acronym": "gene"})
G.to_csv(f"{A}/bitarello_candidate_genes.tsv", sep="\t", index=False)

print(f"S5 windows: {len(W):,} rows  ({(W.set=='OUTLIERS').sum():,} outlier, {(W.set=='SIGNIF').sum():,} signif)")
print(f"   unique OUTLIER windows: {W[W.set=='OUTLIERS'][['chrom','start','end']].drop_duplicates().shape[0]:,}")
print(f"   unique SIGNIF  windows: {W[W.set=='SIGNIF'][['chrom','start','end']].drop_duplicates().shape[0]:,}")
print(f"S4 candidate genes: {len(G)} rows, {G.gene.nunique()} unique genes "
      f"({G[G.set=='outliers'].gene.nunique()} outlier, {G[G.set=='significant'].gene.nunique()} significant)")
print(f"wrote bitarello_windows.hg19.tsv, bitarello_candidate_genes.tsv")
