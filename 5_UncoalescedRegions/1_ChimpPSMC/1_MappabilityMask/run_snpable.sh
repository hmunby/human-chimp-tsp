#!/usr/bin/env bash
# Build the panTro6 mappability mask with SNPable (https://lh3lh3.users.sourceforge.net/snpable.shtml).
#
# Every overlapping k-mer of the reference is realigned to the reference; a position is callable
# if a high enough fraction of the k-mers covering it map back uniquely and without mismatch.
# The fraction is the stringency r. Two masks are used downstream by ../2_Multihetsep/:
#   r = 0.50  the main analysis
#   r = 0.75  the high-stringency sensitivity analysis
#
# Usage: run_snpable.sh <reference.fa> <outdir>
# Requires SNPable (splitfa, gen_raw_mask.pl, gen_mask) and bwa on PATH.
#
#SBATCH --account=palab
#SBATCH --job-name=snpable
#SBATCH -c 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem-per-cpu=4GB
#SBATCH --time 120:00:00
#SBATCH -o logs/%A.out
#SBATCH -e logs/%A.err
set -euo pipefail

GENOME=${1:?usage: run_snpable.sh <reference.fa> <outdir>}
OUTDIR=${2:?usage: run_snpable.sh <reference.fa> <outdir>}
PREFIX=PanTro6_final_mappability
K=150
THREADS=8

mkdir -p "$OUTDIR"
cd "$OUTDIR"

# 1. extract every overlapping K-mer as a read
echo "extracting overlapping ${K}-mers"
splitfa "$GENOME" $K | split -l 20000000
cat x* > ${PREFIX}_split.${K}
rm -f x*

# 2. align them back to the reference
if [ ! -f "${GENOME}.bwt" ]; then
    echo "indexing $GENOME"
    bwa index "$GENOME"
fi
echo "aligning ${K}-mers"
bwa aln -t $THREADS -R 1000000 -O 3 -E 3 "$GENOME" ${PREFIX}_split.${K} > ${PREFIX}_split.${K}.sai
bwa samse -f ${PREFIX}_split.${K}.sam "$GENOME" ${PREFIX}_split.${K}.sai ${PREFIX}_split.${K}

# 3. raw mask, then one final mask per stringency
echo "generating raw mask"
gen_raw_mask.pl ${PREFIX}_split.${K}.sam > ${PREFIX}_rawMask.${K}.fa
for spec in "0.5 50" "0.75 75"; do
    set -- $spec
    r=$1; tag=$2
    out=${PREFIX}_mask.${K}.${tag}.fa
    echo "generating final mask at r=${r} -> ${out}"
    gen_mask -l $K -r "$r" ${PREFIX}_rawMask.${K}.fa > "$out"
done

echo "done; convert the mask FASTAs to BED with generate_bed.py"
