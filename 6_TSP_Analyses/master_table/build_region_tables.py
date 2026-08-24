#!/usr/bin/env python3
"""Table 1 and its supplementary all-regions equivalent, for any TSP set.

Two tables with identical columns, differing only in which regions they include:

  Table 1              regions holding MORE THAN ONE TSP -- the robust core of the set
  Supplementary        EVERY reported region, singletons included

Both are written as markdown (to paste into the manuscript) and as TSV and CSV (to paste
into Excel / Sheets). Nothing here is recomputed: every column is read from the output of
the stage that owns it, so the tables cannot drift from the analyses they summarise.

| column | from |
|---|---|
| Region (hg19), TSPs, Physical span | `../../4_SINGER.../2_ClusteringLD/tsp_clusters_10kb.<set>.tsv` |
| Gene, Assignment | `../../4_SINGER.../3_GeneAssignment/region_gene_assignments.<set>.tsv` |
| Genetic span | `../../4_SINGER.../2_ClusteringLD/tsp_region_genetic_span.tsv` |
| Variant type | snpEff `Annotation` on the TSPs, from `tsps.<set>.hg19.txt` |
| BMAP decile | region mean B, binned on genome-wide 1 kb-window deciles |
| Literature | Leffler 2013, DeGiorgio 2014, Rasmussen 2014, Bitarello 2018, Rocha 2026 |
| GWAS (L2G effector) | `l2g_genes` from the gene assignment, plus `tsp_all_gwas_hits.long.tsv` |

Literature matching follows the manuscript's stated criteria, which differ by study:
  Leffler / Rasmussen / Bitarello  coordinate overlap with a reported interval or window
  DeGiorgio / Rocha               the region's TSPs fall in the transcribed span of a
                                  listed gene -- these studies rank genes, not intervals.
                                  Rocha's list is the published supplementary 6-Mya table
                                  restricted to genes old in BOTH human and chimpanzee
                                  (population-resolved); see rocha_genes().
                                  Applied independently of the gene assigned here, so a
                                  region can match on a gene it is not assigned to (this
                                  is why SH3YL1 carries a DeGiorgio dagger: its TSPs lie
                                  in SNTG2).

Usage:
    python build_region_tables.py --set pooled
    python build_region_tables.py --set strict --out-dir tables_strict
"""
import argparse
import os
import re

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "../.."))
CL = f"{REPO}/4_SINGER_GenealogicalReconstruction_TSP/4_TSPRegions/2_ClusteringLD"
GA = f"{REPO}/4_SINGER_GenealogicalReconstruction_TSP/4_TSPRegions/3_GeneAssignment"
LIT = f"{REPO}/6_TSP_Analyses/literature"

# Genome-wide 1 kb-window B-value decile edges (upper bound of deciles 1-9).
# 1 = strongest background selection, 10 = weakest.
BMAP_DECILE_EDGES = [607, 695, 757, 803, 839, 869, 895, 922, 949]

# GENCODE transcribed spans, for the gene-based literature criteria. Built by
# build_gene_spans.py, which excludes read-through transcripts -- GENCODE's raw `gene`
# records span every transcript, so HBE1 reads as 237 kb rather than 1.6 kb and swallows ten
# olfactory receptor genes. That inflation alone produced a spurious Rocha citation on the
# OR51E2 region. See build_gene_spans.py for why a length cutoff cannot substitute.
GENCODE_BED = f"{HERE}/resources/" \
              "gencode.v49.genes.protein_coding.autosomes.no_readthrough.bed"

# Regions that within_region_ld.py shows are NOT one linked block are reported as their
# separate LD groups, because one span across two unlinked blocks describes a distance no
# observed haplotype spans. Keyed by region id -> [(start, end), ...].
LD_SUBGROUPS = {
    "4_57918296_57919705": [(57918296, 57918492), (57919221, 57919705)],   # IGFBP7
}

# Curated literature calls that gene-body containment cannot make correctly.
#
# The glycophorin locus is a cluster of paralogues (FREM3, GYPE, GYPB, GYPA) in which the
# trans-species variation is intergenic: our TSPs sit 34 kb from FREM3 and 130 kb from
# GYPE, inside no gene body. Rocha et al. report this locus under those gene names, and it
# is unambiguously the same region, so containment would drop a true match. Padding the
# criterion far enough to catch it would add four spurious DeGiorgio regions, so the
# exception is made explicitly here instead.
#
#   region id -> {study: reason}
LITERATURE_OVERRIDES = {
    "4_144655489_144662054": {
        "Rocha et al. 2026":
            "glycophorin locus (FREM3/GYPE/GYPB/GYPA); TSPs are intergenic within the "
            "cluster, 34 kb from FREM3",
    },
}

COLUMNS = ["Region (hg19)", "Gene", "Assignment", "TSPs", "Physical span (bp)",
           "Genetic span (x10-3 cM)", "Variant type", "BMAP decile", "Literature",
           "GWAS (L2G effector)"]


# Study-technical annotations that identify the STUDY, not the trait. Stripping these
# collapses "Calcium levels" and "Calcium levels (UKB data field 30680)" -- the same trait
# reported twice -- which otherwise crowd out genuinely distinct traits in a truncated
# cell. Parentheticals that DISTINGUISH traits (e.g. "(femur neck)" vs "(femur shaft)")
# are deliberately left alone.
_TECHNICAL = re.compile(
    r"\s*\((?:"
    r"UKB data field[^)]*|PheCode[^)]*|MTAG|PGS-adjusted|baseline|2df test|"
    r"rank-based inverse normal transformation|"
    r"[A-Za-z0-9]{2,8},[^)]*(?:inv-norm|inverse normal)[^)]*|"
    r"[^)]*inv-norm(?:al)? transformed|"
    r"[A-Z][A-Za-z0-9]{1,6}"          # a bare acronym for the trait itself, e.g. (BMI), (ED)
    r")\)", re.I)


def normalise_trait(t):
    t = _TECHNICAL.sub("", str(t))
    return re.sub(r"\s{2,}", " ", t).strip(" ,;")


def format_traits(traits, keep=2):
    """Distinct traits, most-supported first, truncated for the table cell."""
    counts = {}
    for t in traits:
        n = normalise_trait(t)
        if n:
            counts[n] = counts.get(n, 0) + 1
    order = sorted(counts, key=lambda k: (-counts[k], k.lower()))
    if not order:
        return "—"
    shown = "; ".join(order[:keep])
    return shown + (f" (+{len(order) - keep} more)" if len(order) > keep else "")


def load_literature():
    lit = {}
    lef = pd.read_csv(f"{LIT}/leffler_regions/leffler_125_regions.hg19.tsv", sep="\t")
    lef = lef[lef["resequenced"].isin(["yes", "no"])]
    lit["Leffler et al. 2013"] = ("interval", lef.rename(
        columns={"start_hg19": "start", "end_hg19": "end"})[["chrom", "start", "end"]])

    ras = pd.read_csv(f"{LIT}/rasmussen_2014/rasmussen_2014_top20_highTMRCA_regions.hg19.tsv",
                      sep="\t", comment="#")
    ras["chrom"] = ras["chrom"].astype(str).str.replace("chr", "", regex=False)
    lit["Rasmussen et al. 2014"] = ("interval", ras.rename(
        columns={"start_hg19": "start", "end_hg19": "end"})[["chrom", "start", "end"]])

    bit = pd.read_csv(f"{LIT}/bitarello_2018/bitarello_windows.hg19.tsv", sep="\t")
    bit = bit[bit["set"] == "OUTLIERS"][["chrom", "start", "end"]].drop_duplicates()
    lit["Bitarello et al. 2018"] = ("interval", bit)

    deg = {l.strip().upper() for l in open(f"{LIT}/degiorgio_2014/degiorgio_top100_genes.txt")
           if l.strip()}
    lit["DeGiorgio et al. 2014"] = ("gene", deg)
    lit["Rocha et al. 2026"] = ("gene", rocha_genes())
    return lit


# Rocha's 6-Mya outlier genes, taken from their published supplementary tables. The
# criterion is the manuscript's: genes reported as deeply coalescing in BOTH humans and
# chimpanzees. Rocha resolve this at two levels and they disagree, because pooling the chimp
# populations dilutes a signal private to one or two of them:
#
#   species-level      8 genes   HBE1, SNTG2 among our regions
#   population-resolved 29 genes  adds GPR158 and IGFBP7, both `human_only` once the chimp
#                                 populations are pooled but chimp+human when resolved
#
# Population-resolved is used, matching the results text. `all_three` (bonobo + chimp +
# human) is included: those genes are shared with chimpanzee too.
ROCHA_LEVEL = "pops"          # "pops" (reported) | "species"
ROCHA_SHARING = ("chimp_human", "all_three")


def rocha_genes(level=None):
    """Rocha genes deeply coalescing in both humans and chimpanzees."""
    f = ("SuppTable_species_pops_overlap_summary_6MYA_NatGenetics.csv"
         if (level or ROCHA_LEVEL) == "pops" else
         "SuppTable_species_overlap_summary_6MYA_NatGenetics.csv")
    d = pd.read_csv(f"{LIT}/rocha_2026/{f}")
    return {g.strip().upper() for g in
            d.loc[d["species_sharing"].isin(ROCHA_SHARING), "genes_renamed"].astype(str)}


def gene_spans():
    """Transcribed span of each protein-coding gene, for the gene-based criteria."""
    if not os.path.exists(GENCODE_BED):
        print(f"  note: {os.path.basename(GENCODE_BED)} absent -- gene-based literature "
              "criteria (DeGiorgio, Rocha) will be skipped")
        return None
    g = pd.read_csv(GENCODE_BED, sep="\t", header=None, names=["chrom", "start", "end", "gene"])
    g["chrom"] = g["chrom"].astype(str).str.replace("chr", "", regex=False)
    g["gene"] = g["gene"].astype(str).str.upper()
    return g


def literature_for(region_row, tsps, lit, genes):
    """Which prior studies report this region, per the manuscript's stated criteria.

    Interval criteria use hg19 (the coordinate space the published intervals are in);
    the gene-body criteria use hg38, matching GENCODE v49.
    """
    hits, notes = [], []
    c, s, e = str(region_row["chrom"]), int(region_row["start"]), int(region_row["end"])
    for name, (kind, ref) in lit.items():
        if kind == "interval":
            r = ref[ref["chrom"].astype(str) == c]
            if len(r) and ((r["start"] <= e) & (r["end"] >= s)).any():
                hits.append(name)
        elif genes is not None:
            # the region's TSPs must fall inside a listed gene's transcribed span
            gsub = genes[(genes["chrom"] == c) & (genes["gene"].isin(ref))]
            if len(gsub):
                # GENCODE is hg38, so match on the TSPs' hg38 coordinates
                pos = (tsps.dropna(subset=["POSITION_HG38"])["POSITION_HG38"]
                       .astype(int).to_numpy())
                inside = gsub[gsub.apply(
                    lambda g: ((pos >= g["start"]) & (pos <= g["end"])).any(), axis=1)]
                if len(inside):
                    hits.append(name)
                    if inside["gene"].iloc[0].upper() != str(region_row["assigned_gene"]).upper():
                        notes.append(f"{name.split()[0]}: TSPs lie in {inside['gene'].iloc[0]}")
    return ("; ".join(hits) if hits else "—"), "; ".join(notes)


def build(setname, out_dir):
    cl = pd.read_csv(f"{CL}/tsp_clusters_10kb.{setname}.tsv", sep="\t")
    tsps = pd.read_csv(f"{CL}/tsps.{setname}.hg19.txt", sep="\t", low_memory=False)
    # hg38 columns are POSITION_HG38 from aggregate_seeds.py but POS_HG38 from some
    # upstream builders; normalise so the GWAS coordinate match works either way
    if "POSITION_HG38" not in tsps.columns and "POS_HG38" in tsps.columns:
        tsps = tsps.rename(columns={"POS_HG38": "POSITION_HG38"})
    ga = pd.read_csv(f"{GA}/region_gene_assignments.{setname}.tsv", sep="\t").set_index("region")
    span = pd.read_csv(f"{CL}/tsp_region_genetic_span.tsv", sep="\t").set_index("region")
    # region_genetic_span.py also emits a comparison table in which any region that is not
    # one LD block appears as its separate blocks; that is where the per-block cM spans come
    # from, so a split region reports a genetic span per block just as it does a physical one
    block_cm = {}
    cmp_f = f"{CL}/tsp_region_genetic_span_comparison.tsv"
    if os.path.exists(cmp_f):
        cmp_df = pd.read_csv(cmp_f, sep="\t")
        # explicit column names: itertuples mangles "Region (hg19)" and
        # "Span_1e-3cM_HapMapII" into positional attributes, which silently returned the
        # default and left split regions showing a single whole-region span
        loc_col = "Region (hg19)"
        cm_col = next(c for c in cmp_df.columns if "HapMapII" in c)
        for _, x in cmp_df.iterrows():
            loc = str(x[loc_col])                      # "chr4:57,918,296-57,918,492"
            if ":" not in loc:
                continue
            co, rng = loc.replace("chr", "").split(":")
            a, b = (int(v.replace(",", "")) for v in rng.split("-"))
            block_cm[(co, a, b)] = float(x[cm_col])
    gwas_f = f"{GA}/tsp_all_gwas_hits.long.tsv"
    gwas_by_pos = {}
    if os.path.exists(gwas_f):
        gw_df = pd.read_csv(gwas_f, sep="\t")
        for _, x in gw_df.iterrows():
            vid = str(x.get("variant_id_hg38", ""))
            parts = vid.split("_")
            if len(parts) >= 2 and str(x.get("gwas_trait", "")).strip():
                gwas_by_pos.setdefault((parts[0], int(parts[1])), []).append(x["gwas_trait"])
    lit, genes = load_literature(), gene_spans()

    rows, footnotes = [], []
    for r in cl.sort_values("n_tsps", ascending=False).itertuples():
        t = tsps[tsps["region"] == r.region]
        g = ga.loc[r.region] if r.region in ga.index else None

        # physical / genetic span, split where the region is not one LD block
        if r.region in LD_SUBGROUPS:
            blks = LD_SUBGROUPS[r.region]
            phys = " / ".join(str(b - a) for a, b in blks)
            per = [block_cm.get((str(r.chrom), a, b)) for a, b in blks]
            if all(v is not None and pd.notna(v) for v in per):
                gen = " / ".join(f"{v:.3f}" for v in per) + "*"
            else:   # comparison table absent: fall back to the whole-region span
                sp = span.loc[r.region, "span_cM"] * 1e3 if r.region in span.index else np.nan
                gen = f"{sp:.3f}*"
            footnotes.append(f"*{g['assigned_gene'] if g is not None else r.region} comprises "
                             f"{len(blks)} distinct linkage-disequilibrium groups rather than "
                             "one contiguous haplotype; physical and genetic spans are given "
                             "for each group.")
        else:
            phys = str(int(r.span_bp))
            gen = (f"{span.loc[r.region, 'span_cM'] * 1e3:.3f}"
                   if r.region in span.index else "—")

        ann = t["Annotation"].astype(str).str.replace("_variant", "", regex=False)
        vt = ", ".join(sorted(set(ann.str.replace("_region", "", regex=False)))) or "—"

        b = pd.to_numeric(t["BMAP"], errors="coerce").mean()
        dec = int(np.digitize([b], BMAP_DECILE_EDGES, right=True)[0] + 1) if pd.notna(b) else 0

        extra = LITERATURE_OVERRIDES.get(r.region, {})
        litstr, note = literature_for(
            {"chrom": r.chrom, "start": r.start, "end": r.end,
             "assigned_gene": g["assigned_gene"] if g is not None else ""}, t, lit, genes)
        if extra:
            have = [] if litstr == "—" else litstr.split("; ")
            for st, why in extra.items():
                if st not in have:
                    have.append(st)
                    footnotes.append(f"‡{st.split()[0]}: {why}")
            litstr = "; ".join(sorted(have))
        if note:
            footnotes.append(f"†{note}")

        # GWAS TRAITS reported at this region's TSPs. No gene fallback: this column is
        # traits, and l2g_genes (the assignment evidence) belongs in the Gene column.
        tr = []
        for x in t.dropna(subset=["POSITION_HG38"]).itertuples():
            tr += gwas_by_pos.get((str(x.CHROM_HG38), int(x.POSITION_HG38)), [])
        gw = format_traits(tr)

        rows.append({
            "Region (hg19)": f"chr{r.chrom}:{r.start:,}-{r.end:,}",
            "Gene": g["assigned_gene"] if g is not None else "—",
            "Assignment": g["assignment_tier"] if g is not None else "—",
            "TSPs": int(r.n_tsps),
            "Physical span (bp)": phys,
            "Genetic span (x10-3 cM)": gen,
            "Variant type": vt,
            "BMAP decile": dec,
            "Literature": litstr,
            "GWAS (L2G effector)": gw,
        })

    df = pd.DataFrame(rows)[COLUMNS]
    os.makedirs(out_dir, exist_ok=True)
    made = []
    for tag, sub in (("Table1_multiTSP", df[df["TSPs"] > 1]),
                     ("SuppTable_all_regions", df)):
        base = f"{out_dir}/{tag}.{setname}"
        sub.to_csv(f"{base}.tsv", sep="\t", index=False)
        sub.to_csv(f"{base}.csv", index=False)
        with open(f"{base}.md", "w") as fh:
            title = ("Table 1 — Regions with multiple trans-species polymorphisms"
                     if tag.startswith("Table1") else
                     "Supplementary Table — All trans-species polymorphic regions")
            fh.write(f"# {title}\n\n({setname} TSP set: {len(df)} regions, "
                     f"{int(df['TSPs'].sum())} TSPs; this table shows {len(sub)})\n\n")
            fh.write(sub.to_markdown(index=False) + "\n")
            if footnotes:
                fh.write("\n" + "\n".join(f"- {x}" for x in dict.fromkeys(footnotes)) + "\n")
        made.append((tag, len(sub)))
        print(f"  {tag:<22} {len(sub):>3} rows -> {os.path.basename(base)}.{{md,tsv,csv}}")
    return df


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--set", default="pooled", help="TSP set tag (pooled, strict, ne20k, ne100k)")
    ap.add_argument("--out-dir", default=HERE)
    a = ap.parse_args()
    print(f"building region tables for the '{a.set}' set")
    build(a.set, a.out_dir)


if __name__ == "__main__":
    main()
