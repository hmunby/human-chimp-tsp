#!/usr/bin/env bash
# liftOver the minimal-decode beds from panTro6 to hg38, then clean them.
#
# Runs after make_minimal_decode.py and produces the per-chromosome hg38 beds that
# ../../2_UncoalescedSegments/ consumes:
#   1. liftOver each panTro6 bed to hg38 (samples in parallel)
#   2. concatenate and sort per sample
#   3. keep only segments that survived the lift at close to their original 1 kb length
#
# Usage: lift_clean_minimal_decode.sh [minimal_decode_dir] [chain]
set -euo pipefail

DIR=${1:-minimal_decode/5.5Mya}
CHAIN=${2:-resources/panTro6ToHg38.over.chain}
CLEAN="$(dirname "$0")/clean_lifted_minimal_decode.py"

SAMPLES="CEN_SAMEA4374772 CEN_SAMEA4374778 CEN_SAMEA4374785 EAS_2003 EAS_SAMEA4374790 EAS_SAMEA4374797 WES_1059 WES_SAMEA2421542 WES_SAMEA5204228"
CHROMS="1 2A 2B 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"

mkdir -p "$DIR/lift_tmp" "$DIR/cleaned_lift_hg38"

lift_sample () {
  local popsam=$1
  for chrom in $CHROMS; do
    liftOver "$DIR/${popsam}_chr${chrom}_minimal_decode.bed" "$CHAIN" \
      "$DIR/lift_tmp/${popsam}_chr${chrom}_minimal_decode_hg38.bed" \
      "$DIR/lift_tmp/${popsam}_chr${chrom}_minimal_decode_hg38.bed.unmapped" 2>/dev/null
  done
  cat "$DIR"/lift_tmp/${popsam}_chr*_minimal_decode_hg38.bed \
    | sort -k1,1V -k2,2n -k3,3n > "$DIR/${popsam}_minimal_decode_hg38.bed"
  echo "  lifted and combined $popsam"
}
export -f lift_sample
export DIR CHAIN CHROMS

printf '%s\n' $SAMPLES | xargs -P 9 -I{} bash -c 'lift_sample "$@"' _ {}

for popsam in $SAMPLES; do
  python3 "$CLEAN" "$DIR/${popsam}_minimal_decode_hg38.bed" "$DIR/cleaned_lift_hg38" "$popsam"
done

echo "done; cleaned beds in $DIR/cleaned_lift_hg38/"
