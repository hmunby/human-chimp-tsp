# Branch-midpoint age-estimator bias — Supplementary Figure 3

Allele ages are summarised by the midpoint of the branch the mutation is inferred to have arisen
on. The expected age of a neutral mutation given its branch *is* the branch midpoint, but once you
condition on mutations that are genuinely old that stops being true: old mutations fall in the
older half of the deepest branches, and the deepest branches are the longest. 

This quantifies that bias by simulation: neutral coalescent genealogies with the true age of every
mutation known, plotted against what the midpoint estimator would have reported.

## Simulation

Constant-size neutral coalescent under msprime: 2,000 diploid samples (4,000 haplotypes),
Nₑ = 10,000, μ = 1e-8 per site per generation, 1 Mb segments without recombination. Mutations on
the root branch and below 1% derived-allele frequency are excluded. 100 runs of 10⁶ mutations each
(about 5 × 10⁴ replicates in total) give roughly 10⁸ mutations.

The estimator is unbiased for young mutations and bends below the identity line for mutations older
than about 3–4 Nₑ generations.

## Flow (rules, in order)
1. **simulate_neutral** (`scripts/simulate_neutral.py`) — one run per seed, recording each
   mutation's true age and its branch midpoint → `results/neutral_sim_{run}.pkl`.
2. **combine_sims** — concatenate the 100 runs.
3. **compare_ages** (`scripts/compare_ages.py`) — the diagnostic comparison.

Then `scripts/plot_suppfig_midpoint_bias.py` draws the supplementary figure: a random subsample of
10,000 mutations, the identity line, and a loess of the midpoint estimate on the true age computed
over all simulated mutations. `scripts/midpoint_curve_cache.npz` caches that loess so the figure
can be redrawn without the 800 MB combined pickle.

```
snakemake --profile .                      # simulate (long; 100 jobs)
python scripts/plot_suppfig_midpoint_bias.py
```

## Software
Python 3 with msprime, numpy, pandas, matplotlib, statsmodels (loess); Snakemake 7.
