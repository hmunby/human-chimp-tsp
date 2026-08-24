# 1_ChimpPSMC / 3_Cobraa — cobraa fit and posterior decode

Fits cobraa to each individual's multihetsep files and then decodes the posterior along each
chromosome.

## Rules
- **run_cobraa** — one fit per individual across all 23 autosomes, from the r = 0.50 mappability
  mask. This is the main analysis.
- **run_cobraa_high_stringency** — the same fit from the r = 0.75 mask, into `results_75/`, as a
  mappability sensitivity check.
- **decode_cobraa** — re-runs cobraa for a single iteration per chromosome with θ, ρ and λ_A held
  at their fitted values (`-decode -decode_downsample 10`), emitting the posterior over the 64
  coalescent-time intervals every 1 kb. This is the input to `../4_MinimalDecode/`.

## Parameters (published set)

| | |
|---|---|
| time intervals | D = 64 |
| bin size | b = 100 bp |
| spread₁ / spread₂ | 0.075 / 100 |
| initial μ/ρ | 1.5 |
| EM iterations | 30 |
| θ | fixed at 0.001 |
| decode downsampling | 10 (posterior every 1 kb) |

## Inputs to set
- `MHS_DIR` — `../2_Multihetsep/` output root (pre-filled).
- `COBRAA` — path to `cobraa.py` ([cobraa](https://github.com/TrevorCousins/cobraa)) 