#!/usr/bin/env bash
# Segmental-duplication track, step 3: merge overlapping SD intervals into the final mask BED.
# The no-chr version (t2t_pantro6_seg_dups_resolved_no_chr.merged.bed) is the file consumed by the
# filtering Snakemake (rule `segdup`: bcftools view -T ^<bed>)
# Tool: bedtools v2.29.1.
set -euo pipefail

BED_DIR="processed_bed"

for prefix in t2t_pantro6_seg_dups_resolved t2t_pantro6_seg_dups_resolved_no_chr; do
    bedtools merge -i "${BED_DIR}/${prefix}.bed" > "${BED_DIR}/${prefix}.merged.bed"
done

# total masked bp
awk '{sum += $3 - $2} END {print "merged SD bp:", sum}' \
    "${BED_DIR}/t2t_pantro6_seg_dups_resolved.merged.bed"
