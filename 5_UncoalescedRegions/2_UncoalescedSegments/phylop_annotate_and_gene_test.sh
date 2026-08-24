#!/usr/bin/env bash
# Uncoalesced regions: phyloP annotation, the "unconstrained" filter, and the
# closest-gene test set. Sits between build_chimp_and_intersect.py (which writes the chimp
# segments with a phyloP PLACEHOLDER) and the figure and enrichment scripts.
#
# Two INDEPENDENT branches. They do not meet, and neither is derived from the other:
#
#   A. phyloP -> unconstrained -> reciprocal unconstrained intersect
#      Restricts both species to the least-conserved 90% of the genome, so that reciprocal
#      uncoalescence is measured away from constrained sequence.   -> Figure 6
#
#   B. merge the exact reciprocal intersect -> nearest gene
#      Operates on the UNFILTERED exact intersect, not on branch A's output.   -> enrichment
#
# intermediates. See the notes at the bottom for the three things that are easy to get wrong.
set -euo pipefail

# ---- inputs (edit) ------------------------------------------------------------------------------
OUT="expected_uncoal"                     # output tree (build_chimp_and_intersect.py)

# chimp genome-wide segments: chrom,start,end,p0,phylop(PLACEHOLDER),bmap
CHIMP_GW="$OUT/chimp/expected_uncoal_chimp_5755625.genome_wide.bmap.bed"

# human x chimp exact overlap, already written by build_chimp_and_intersect.py
EXACT="$OUT/intersect/expected_uncoal_intersect.whole_genome.exact.bed"

# human "unconstrained" segments
HUMAN_UNCONSTRAINED="expected_uncoal/human/phylop/expected_uncoal_human_5589889.genome_wide.phylop100way.unconstrained_cutoff_0.331998.bed"

BIGWIG="resources/hg38.phyloP100way.bw"         # UCSC hg38 phyloP 100-way bigWig (~9.9 GB)
GENCODE="resources/gencode.v49.pc_genes.bed"    # GENCODE v49 protein-coding gene bodies + names
CUTOFF=0.331998                                 # 90th percentile of genome-wide 1 kb mean phyloP
BIGWIGAVG="bigWigAverageOverBed"
BEDTOOLS="bedtools"
# ------------------------------------------------------------------------------------------------
mkdir -p "$OUT/tmp" "$OUT/chimp" "$OUT/intersect"

CHIMP_PHYLOP="$OUT/chimp/expected_uncoal_chimp_5755625.genome_wide.phylop100way.bmap.bed"
CHIMP_UNCONSTRAINED="$OUT/chimp/expected_uncoal_chimp_5755625.genome_wide.phylop100way.unconstrained_cutoff_${CUTOFF}.bed"
RECIP="$OUT/intersect/expected_uncoal_intersect.unconstrained_cutoff_${CUTOFF}.bed"
MERGED="$OUT/intersect/expected_uncoal_intersect.whole_genome.exact.merged.bed"
CLOSEST="$OUT/intersect/expected_uncoal_intersect.whole_genome.exact.merged.closest_genes.bed"

### BRANCH A -- phyloP, unconstrained filter, reciprocal intersect (Figure 6) ###

# A1. mean phyloP per chimp segment. The name column is the row number so -bedOut can paste the
#     score straight back in row order; the BMAP column is re-attached from the source file.
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,NR}' "$CHIMP_GW" > "$OUT/tmp/chimp_id.bed"
"$BIGWIGAVG" "$BIGWIG" "$OUT/tmp/chimp_id.bed" "$OUT/tmp/chimp.tab" \
    -bedOut="$OUT/tmp/chimp_phylop.bed"
paste "$OUT/tmp/chimp_phylop.bed" <(cut -f6 "$CHIMP_GW") > "$CHIMP_PHYLOP"
# -> chrom, start, end, row_id, mean phyloP, bmap

# A2. chimp "unconstrained" = mean phyloP strictly below the cutoff
awk -F'\t' -v c="$CUTOFF" '$5 < c' "$CHIMP_PHYLOP" > "$CHIMP_UNCONSTRAINED"

# A3. reciprocal unconstrained intersect. -wao keeps every human segment, with zero-fill where no
#     chimp segment overlaps, and appends the overlap length; the figure needs that full 13-column
#     table, not just the overlapping rows.
"$BEDTOOLS" intersect -a "$HUMAN_UNCONSTRAINED" -b "$CHIMP_UNCONSTRAINED" -wao -f 0.0001 -r \
    > "$RECIP"
# -> 6 human cols, 6 chimp cols, overlap_len

### BRANCH B -- merged exact intersect, nearest gene (enrichment) ###

# B1. length-weighted merge of the EXACT intersect (see note 2)
python3 merge_bed_weighted.py "$EXACT" -o "$MERGED"

# B2. nearest protein-coding gene, with distance (see note 3)
sort -k1,1 -k2,2n "$MERGED" > "$OUT/tmp/merged.sorted.bed"
"$BEDTOOLS" closest -a "$OUT/tmp/merged.sorted.bed" -d \
    -b <(sort -k1,1 -k2,2n "$GENCODE") > "$CLOSEST"
# -> 6 segment cols, 4 gene cols, distance

echo "done:"
echo "  Figure 6 input   -> $RECIP"
echo "  enrichment input -> $CLOSEST"