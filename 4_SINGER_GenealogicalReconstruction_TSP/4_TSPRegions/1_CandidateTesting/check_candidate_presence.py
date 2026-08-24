#!/usr/bin/env python3
"""Find regions whose candidate SNP was lost in chimpanzee recalling.

SINGER's downsample step anchors on the candidate SNP, so a region whose candidate did not
survive recalling and filtering fails with "No shared sites between sites and shared_pos files"
and blocks aggregation. This is the pre-flight check: run it after staging the merged region
VCFs, before running SINGER.

A region survives if at least one position in its merged VCF matches a candidate SNP. The check
is done two ways:

  CHROM+POS  the correct test, and what this script excludes on.
  POS only   what downsample itself matches on. A region can therefore be rescued by a
             position that happens to coincide with a candidate on a DIFFERENT chromosome.
             Those are reported separately as "masked-collision": downsample will not crash,
             but the region has no real candidate and should still be excluded.

Run from this directory. Prints the regions to exclude; with --write appends the newly lost ones
to the SINGER workdir's misc_files/excluded_regions.txt, de-duplicated.
"""
import glob
import os
import subprocess
import sys
from collections import defaultdict

# ---- inputs (edit) -------------------------------------------------------------------------
WORK = "../../3_SINGER"   # the SINGER working directory
BCFTOOLS = "bcftools"
# --------------------------------------------------------------------------------------------

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, WORK)
MERGED = os.path.join(WORK, "merged_regions")
CAND = os.path.join(WORK, "misc_files/merged_candidate_snps.chrom_pos.txt")
EXCL = os.path.join(WORK, "misc_files/excluded_regions.txt")

write = "--write" in sys.argv

candidates_by_chrom = defaultdict(set)
candidates_any_chrom = set()
with open(CAND) as f:
    for line in f:
        p = line.split()
        if len(p) < 2 or p[0] == "X":
            continue
        candidates_by_chrom[p[0]].add(int(p[1]))
        candidates_any_chrom.add(int(p[1]))

vcfs = sorted(glob.glob(os.path.join(MERGED, "*.vcf.gz")))
lost = []       # no candidate for this chromosome present
crashes = []    # not even a position-only match, so downsample will crash
masked = []     # rescued only by a cross-chromosome position collision
ok = 0

for v in vcfs:
    region = os.path.basename(v)[:-7]
    chrom = region.split("_")[0]
    out = subprocess.run([BCFTOOLS, "query", "-f", "%POS\n", v],
                         capture_output=True, text=True)
    present = {int(x) for x in out.stdout.split()}
    if present & candidates_by_chrom.get(chrom, set()):
        ok += 1
        continue
    lost.append(region)
    (masked if present & candidates_any_chrom else crashes).append(region)

print(f"staged regions checked: {len(vcfs)}")
print(f"  candidate present (CHROM+POS): {ok}")
print(f"  candidate LOST (CHROM+POS):    {len(lost)}")
print(f"    of which downsample WILL CRASH (no position-only match): {len(crashes)}")
print(f"    of which masked by a cross-chromosome collision:         {len(masked)}")

existing = set()
if os.path.exists(EXCL):
    existing = {line.strip() for line in open(EXCL) if line.strip()}
new_to_exclude = [r for r in lost if r not in existing]

print(f"\ncurrently in excluded_regions.txt: {len(existing)}")
print(f"newly lost (not yet excluded):     {len(new_to_exclude)}")
for r in sorted(new_to_exclude):
    print(f"    {r}\t[{'CRASHES' if r in crashes else 'masked-collision'}]")

# regions excluded earlier whose candidate has since come back in a fresh recall
staged = {os.path.basename(v)[:-7] for v in vcfs}
recovered = [r for r in existing if r in staged and r not in lost]
if recovered:
    print(f"\nNOTE: {len(recovered)} excluded region(s) are staged and now have their candidate "
          "present (recovered in a fresh recall):")
    for r in sorted(recovered):
        print(f"    {r}")

if write and new_to_exclude:
    with open(EXCL, "a") as f:
        for r in sorted(new_to_exclude):
            f.write(r + "\n")
    print(f"\nappended {len(new_to_exclude)} region(s) to {EXCL}")
elif write:
    print("\nnothing to append")
