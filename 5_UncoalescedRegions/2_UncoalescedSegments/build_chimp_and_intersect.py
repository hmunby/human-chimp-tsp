#!/usr/bin/env python3
"""Build the chimp genome-wide expected-uncoal file and the human-chimp
intersections, writing the results into expected_uncoal/
tree. 

Chimp segments: sum P(uncoal) across the 9 samples per 1 kb
segment, keep sum>1; merge (bedtools merge -d -1); annotate BMAP (round start to nearest 1000,
look up the decile-bed raw B value). phylop/p0 columns are placeholders (no downstream figure
uses chimp phylop or p0; Fig 6 uses chimp BMAP, Fig 2/4 use human cols + positions).

Outputs (chr order matches the pipeline's genome-wide concat: chr1,chr10,chr11,...):
  expected_uncoal/chimp/expected_uncoal_chimp_5755625.genome_wide.bmap.bed
  expected_uncoal/intersect/expected_uncoal_intersect.whole_genome.exact.bed   (Fig 2)
  expected_uncoal/intersect/expected_uncoal_intersect.whole_genome.bed         (Fig 4, -wao)
"""
import subprocess, os, tempfile, pandas as pd

AP = "resources"
BEDTOOLS = "bedtools"
CHIMP_DIR = "resources/cleaned_lift_hg38"
HUMAN_GW = f"{AP}/expected_uncoal/human/phylop/expected_uncoal_human_5589889.genome_wide.phylop100way.bmap.bed"
BMAP_DECILE = "resources/whole_genome.1kb.bmap_deciles.bed"
OUT = f"{AP}/expected_uncoal"
SAMPLES = ['EAS_2003','EAS_SAMEA4374790','EAS_SAMEA4374797','WES_1059','WES_SAMEA5204228',
           'WES_SAMEA2421542','CEN_SAMEA4374772','CEN_SAMEA4374778','CEN_SAMEA4374785']
# genome-wide concat order used by the pipeline (string sort of chr1..chr22)
CHR_ORDER = sorted(range(1, 23), key=lambda c: str(c))

os.makedirs(f"{OUT}/chimp", exist_ok=True)
os.makedirs(f"{OUT}/intersect", exist_ok=True)

# raw BMAP lookup: rounded-start -> raw B (per chrom)
bd = pd.read_csv(BMAP_DECILE, sep="\t", header=None, names=["c","s","e","bmap","dec"])
bmap_by_pos = {}
for c, s, b in zip(bd["c"], bd["s"], bd["bmap"]):
    bmap_by_pos[(c, s)] = b

tmp = tempfile.mkdtemp()
gw = open(f"{OUT}/chimp/expected_uncoal_chimp_5755625.genome_wide.bmap.bed", "w")
for chrom in CHR_ORDER:
    acc = {}
    for s in SAMPLES:
        f = f"{CHIMP_DIR}/{s}_chr{chrom}_minimal_decode_hg38.bed"
        if not os.path.exists(f):
            continue
        for line in open(f):
            v = line.split("\t")
            key = (int(float(v[1])), int(float(v[2])))
            acc[key] = acc.get(key, 0.0) + float(v[3])
    seg = [(f"chr{chrom}", st, en) for (st, en), p in acc.items() if p > 1]
    if not seg:
        continue
    raw = f"{tmp}/c{chrom}.bed"
    pd.DataFrame(seg, columns=["c","s","e"]).sort_values(["s","e"]).to_csv(raw, sep="\t", header=False, index=False)
    merged = f"{tmp}/c{chrom}.merged.bed"
    with open(merged, "w") as o:
        subprocess.run([BEDTOOLS, "merge", "-d", "-1", "-i", raw], stdout=o, check=True)
    cc = f"chr{chrom}"
    for line in open(merged):
        c, st, en = line.split("\t"); st, en = int(st), int(en)
        rstart = 1000 * round(st / 1000)
        b = bmap_by_pos.get((cc, rstart), "NA")
        gw.write(f"{cc}\t{st}\t{en}\t1.0\t0.0\t{b}\n")
gw.close()
chimp_gw = f"{OUT}/chimp/expected_uncoal_chimp_5755625.genome_wide.bmap.bed"
n = sum(1 for _ in open(chimp_gw))
print(f"chimp genome-wide: {n:,} segments -> {chimp_gw}")

# intersects with the human genome-wide bmap
exact = f"{OUT}/intersect/expected_uncoal_intersect.whole_genome.exact.bed"
wao = f"{OUT}/intersect/expected_uncoal_intersect.whole_genome.bed"
with open(exact, "w") as o:
    subprocess.run([BEDTOOLS, "intersect", "-a", HUMAN_GW, "-b", chimp_gw], stdout=o, check=True)
with open(wao, "w") as o:
    subprocess.run([BEDTOOLS, "intersect", "-a", HUMAN_GW, "-b", chimp_gw, "-wao", "-f", "0.0001", "-r"],
                   stdout=o, check=True)
print(f"exact intersect -> {exact}  ({sum(1 for _ in open(exact)):,} rows)")
print(f"wao   intersect -> {wao}  ({sum(1 for _ in open(wao)):,} rows)")
