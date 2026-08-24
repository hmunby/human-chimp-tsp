#!/usr/bin/env bash
# Per-individual chimp uncoalesced-by-decile 
set -euo pipefail
WD=../../../5_UncoalescedRegions/3_Figures/per_individual
DECILE=resources/whole_genome.1kb.bmap_deciles.bed
CHIMP_DIR=resources/cleaned_lift_hg38
BEDTOOLS=${BEDTOOLS:-bedtools}
MHC_C=chr6; MHC_S=28510120; MHC_E=33480577
CHIMPS="CEN_SAMEA4374772 CEN_SAMEA4374778 CEN_SAMEA4374785 EAS_2003 EAS_SAMEA4374790 EAS_SAMEA4374797 WES_1059 WES_SAMEA2421542 WES_SAMEA5204228"
HUMANS="ESN_HG03515 YRI_NA18488 MSL_HG03212 GWD_HG02568 ACB_HG01882 ASW_NA19625 LWK_NA19017"

accumulate () {
  awk -v C="$MHC_C" -v S="$MHC_S" -v E="$MHC_E" 'BEGIN{OFS="\t"} !($1==C && $3>S && $2<E)' \
  | sort -k1,1 -k2,2n \
  | $BEDTOOLS intersect -a - -b "$DECILE" -wo \
  | awk 'BEGIN{OFS="\t"} {d=$9; ov=$10; num[d]+=$4*ov; den[d]+=ov}
         END{for(i=1;i<=10;i++){f=(den[i]>0)?num[i]/den[i]:0; print i, den[i]+0, num[i]+0, f}}'
}

mkdir -p "$WD/chimp"
for s in $CHIMPS; do
  echo "chimp $s ..."
  cat "$CHIMP_DIR/${s}_chr"{1..22}"_minimal_decode_hg38.bed" | accumulate > "$WD/chimp/${s}.decile_frac.tsv"
done

# assemble the data dirs: chimp built above
for pair in "data:data" "data_sd:data_sd" "data_nosd:data_nosd"; do
  src="${pair%%:*}"; dst="${pair##*:}"
  mkdir -p "$WD/$dst"
  for s in $CHIMPS; do cp "$WD/chimp/${s}.decile_frac.tsv" "$WD/$dst/${s}.decile_frac.tsv"; done
  for s in $HUMANS;  do cp "$WD/$src/${s}.decile_frac.tsv"       "$WD/$dst/${s}.decile_frac.tsv"; done
  echo "assembled $dst ($(ls "$WD/$dst"/*.tsv | wc -l) samples)"
done
echo DONE
