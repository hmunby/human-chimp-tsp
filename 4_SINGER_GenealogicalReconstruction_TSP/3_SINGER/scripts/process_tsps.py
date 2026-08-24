"""Aggregate one run's per-region SINGER metrics into that run's TSP set.

    candidate SNP  ->  support > threshold  ->  TSPs

Steps:
 1. Aggregate every target region's grouped all-shared metrics.
 2. Restrict to the ORIGINAL CANDIDATE SNPs.
 3. Keep those supported as TSP (first- or second-coalescence) in more than
    `tsp_threshold` per cent of the retained posterior ARG samples.
 4. Merge with the shared-SNP annotation table (hg19) and lift to hg38.

The age criterion is NOT applied here. It is what selected the candidates in the first
place (shared SNPs with a YRI midpoint age over 4 Mya), so re-imposing it after the ARG
test would apply the same cut twice.

    Metrics are computed at every shared SNP in a region, not only the candidates, so the
    restriction in step 2 is what makes this set the reported one. It matches
    ../4_TSPRegions/2_ClusteringLD/aggregate_seeds.py, so the single-chain and pooled paths
    agree on what is eligible.

This is the SINGLE-CHAIN call, used directly for the Ne = 20,000 and Ne = 100,000
sensitivity runs. The reported set pools the three Ne = 50,000 chains downstream, in
aggregate_seeds.py, which also does the 10 kb single-linkage clustering into regions.

Run via Snakemake (uses the `snakemake` object).
"""
import glob
import os
import subprocess
import tempfile

import pandas as pd

singer_dir = snakemake.params.singer_dir
cand_f = snakemake.input.candidates
snp_info_f = snakemake.input.snp_info
chain = snakemake.input.chain
thresh = float(snakemake.params.tsp_threshold)
gen = int(snakemake.params.gen_time)

# Current target regions only. singer/ can retain stale region directories from earlier
# windowings that overlap current regions and would double-count loci.
_target = {l.strip() for l in open("misc_files/non_overlapping_region_250kb_ids.txt") if l.strip()}
_excl = set()
if os.path.exists("misc_files/excluded_regions.txt"):
    _excl = {l.split()[0] for l in open("misc_files/excluded_regions.txt")
             if l.strip() and not l.startswith("#")}
_valid = _target - _excl

# ---- 1. aggregate per-region grouped metrics ----------------------------------------------------
dfs = []
for f in glob.glob(os.path.join(singer_dir, "*", "*.all_shared_metrics.grouped.txt")):
    region = os.path.basename(os.path.dirname(f))
    if region not in _valid:
        continue
    # require a currently staged recall VCF, so metrics from a region that has since been
    # re-called but not yet rerun never enter the set
    if not os.path.exists(os.path.join("merged_regions", region + ".vcf.gz")):
        continue
    try:
        d = pd.read_csv(f, sep="\t")
    except Exception:
        continue
    if len(d) == 0:
        continue
    d["region"] = region
    dfs.append(d)

all_df = pd.concat(dfs, ignore_index=True)
all_df = all_df[all_df["chrom"] != "chrom"]                    # the grouped files repeat a header
all_df["chrom"] = all_df["chrom"].astype(int)
all_df["position"] = all_df["position"].astype(int)
all_df = all_df.sort_values(["chrom", "position"])
print(f"regions aggregated: {all_df['region'].nunique()}; shared-SNP rows: {len(all_df)}")

# ---- 2. restrict to the original candidates -----------------------------------------------------
cand = set()
for line in open(cand_f):
    p = line.split()
    if len(p) >= 2 and p[0].isdigit():
        cand.add((int(p[0]), int(p[1])))
print(f"candidate SNPs: {len(cand)}")

key = list(zip(all_df["chrom"], all_df["position"]))
cand_df = all_df[[k in cand for k in key]].copy()
print(f"candidates with metrics in this run: {len(cand_df)} "
      f"({len(cand) - len(cand_df)} not scored; see the stage README on recall dropout)")

# ---- 3. support threshold -----------------------------------------------------------------------
cand_df["tsp_first_or_second_count"] = cand_df["tsp_first_or_second_count"].astype(float)
tsp_df = cand_df[cand_df["tsp_first_or_second_count"] > thresh].copy()
print(f"TSPs (support > {thresh:.0f}%): {len(tsp_df)} in {tsp_df['region'].nunique()} regions")

# ---- 4. annotate and lift -----------------------------------------------------------------------
info = pd.read_csv(snp_info_f, sep="\t", low_memory=False)


def midage(r):
    lo, up = r["LOWER_AGE_YRI"], r["UPPER_AGE_YRI"]
    if pd.isna(lo) or lo == "." or pd.isna(up) or up == ".":
        return float("nan")
    return (float(lo) + float(up)) / 2 * gen


# every candidate with its support, whether or not it passed -- the input to aggregate_seeds.py
allc = pd.merge(info, cand_df, left_on=["CHROM", "POS"], right_on=["chrom", "position"],
                how="right")
allc["MID_AGE_YRI_YEARS"] = allc.apply(midage, axis=1)
allc.to_csv(snakemake.output.all_candidates, sep="\t", index=False)

final = pd.merge(info, tsp_df, left_on=["CHROM", "POS"], right_on=["chrom", "position"],
                 how="right")
final["MID_AGE_YRI_YEARS"] = final.apply(midage, axis=1)

with tempfile.TemporaryDirectory() as td:
    h19, h38, un = (os.path.join(td, x) for x in ("h19.bed", "h38.bed", "un.bed"))
    fb = final[["CHROM", "POS"]].copy()
    fb["start"] = fb["POS"] - 1
    fb["name"] = fb["CHROM"].astype(str) + "_" + fb["POS"].astype(str)
    fb[["CHROM", "start", "POS", "name"]].to_csv(h19, sep="\t", header=False, index=False)
    subprocess.run(["liftOver", h19, chain, h38, un], check=False)
    if os.path.exists(h38) and os.path.getsize(h38) > 0:
        lift = pd.read_csv(h38, sep="\t", header=None, names=["c", "s", "e", "name"])
        lift["CHROM_HG38"] = lift["c"].astype(str).str.replace("chr", "", regex=False)
        lift["POS_HG38"] = lift["e"]
        lmap = lift.set_index("name")[["CHROM_HG38", "POS_HG38"]]
    else:
        lmap = pd.DataFrame(columns=["CHROM_HG38", "POS_HG38"])

final["name"] = final["CHROM"].astype(str) + "_" + final["POS"].astype(str)
final = final.join(lmap, on="name").drop(columns="name")
final.to_csv(snakemake.output.tsps, sep="\t", index=False)

fb = final.dropna(subset=["POS_HG38"]).copy()
if len(fb):
    fb["start"] = fb["POS_HG38"].astype(int) - 1
    fb["end"] = fb["POS_HG38"].astype(int)
    fb[["CHROM_HG38", "start", "end"]].sort_values(["CHROM_HG38", "start"]).to_csv(
        snakemake.output.bed, sep="\t", header=False, index=False)
else:
    open(snakemake.output.bed, "w").close()
print("done")
