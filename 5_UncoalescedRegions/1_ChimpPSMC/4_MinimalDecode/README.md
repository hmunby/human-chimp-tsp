# 1_ChimpPSMC / 4_MinimalDecode — P(uncoalesced at the split), lifted to hg38

Reduces the cobraa posterior grids to one number per 1 kb bin per individual — the probability
that that individual's two haplotypes have **not** coalesced by the human–chimp split — and puts
it in human coordinates.

## Steps
1. `make_minimal_decode.py` — takes the inverse cumulative posterior at the time-grid
   boundary just above 5.5 Mya and writes it as a BED in panTro6 coordinates.
2. `lift_clean_minimal_decode.sh` — liftOver panTro6 → hg38 (samples in parallel), then
   concatenate and sort per sample, then clean.
3. `clean_lifted_minimal_decode.py` — keeps only segments that came through the lift at
   950–1050 bp. A segment far from its original 1 kb straddled a rearrangement or an indel-rich
   region, so its posterior no longer describes one comparable interval.

```
python make_minimal_decode.py
./lift_clean_minimal_decode.sh [minimal_decode_dir] [chain]
```

Output: `minimal_decode/5.5Mya/cleaned_lift_hg38/{popsam}_chr{c}_minimal_decode_hg38.bed`,
read by `../../2_UncoalescedSegments/build_chimp_and_intersect.py`.

## Inputs to set
- `COBRAA_RESULTS` / `PARAMS` in `make_minimal_decode.py` — must match `../3_Cobraa/`.
- the UCSC `panTro6ToHg38.over.chain`, passed as the second argument to the shell script.

## Software
Python 3 (numpy, pandas); UCSC `liftOver`.
