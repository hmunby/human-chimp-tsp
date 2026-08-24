# Per-individual uncoalescence by B-map decile — Supplementary Figure 5

The percentage of each individual's genome that remains uncoalesced at 5.5 Mya, by decile of
background selection. Sixteen individuals: nine chimpanzees (EAS blue, WES purple, CEN red) and
the seven African-ancestry humans of Cousins et al. One line per individual.

## Files
- `build_chimp.sh` — builds the per-sample decile fractions, intersecting the
  cleaned hg38 minimal-decode beds from `../../1_ChimpPSMC/4_MinimalDecode/` with the genome-wide
  1 kb B-map decile bed.
- `plot_per_individual.py` — Takes `DATADIR OUTSUFFIX MODE`, where MODE is `pct` (the
  reported version) or `count`.
- `data/` — the 16 per-sample decile fraction tables, bundled (about 340 B each).

```
python plot_per_individual.py data '' pct
```
## Software
Python 3 (pandas, matplotlib); bedtools.
