#!/usr/bin/env bash
# Download the Roulette per-chromosome mutation-rate VCFs (hg19, v5.2 with TFBS correction) used to
# build the MR annotation. Roulette: Seplyarskiy et al. 2023 (Sunyaev lab).
set -euo pipefail

BASE="http://genetics.bwh.harvard.edu/downloads/Vova/Roulette/hg19/autosomes"
OUT_DIR="resources/roulette/autosomes"
mkdir -p "$OUT_DIR"

for i in $(seq 1 22); do
    f="hg19_chr${i}_rate_v5.2_TFBS_correction_sorted.vcf.gz"
    wget -P "$OUT_DIR" "${BASE}/${f}"
    wget -P "$OUT_DIR" "${BASE}/${f}.tbi"
done
