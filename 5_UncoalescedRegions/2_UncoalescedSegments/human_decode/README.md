# Human PSMC/cobraa decode — the 5.5 Mya column

The human side of the uncoalesced-regions analysis rests on posterior decodings from Cousins et
al. for seven individuals of recent African ancestry. The published decode gives, for every 1 kb
window, the probability of remaining uncoalesced at each of 19 time boundaries older than 2 Mya —
3.4 GB across 154 files. Only one of those boundaries is ever used, so `data/` holds that column
alone and `build_human_uncoalesced.py` rebuilds the human segments from it.

## The column

`PO_5589889` — the probability that an individual's two haplotypes have **not** coalesced by
5,589,889 years, the first boundary past the 5.5 Mya human–chimpanzee split on the human time
grid (µ = 1.25 × 10⁻⁸, generation time 29). The chimpanzee decode uses its own grid and lands on
5,755,625 years instead; see [`../../1_ChimpPSMC/4_MinimalDecode/`](../../1_ChimpPSMC/4_MinimalDecode/).

## `data/` — 7 files, 33 MB

One file per individual: `ACB_HG01882`, `ASW_NA19625`, `ESN_HG03515`, `GWD_HG02568`,
`LWK_NA19017`, `MSL_HG03212`, `YRI_NA18488`. Roughly 2.87 M windows each, all 22 autosomes.

Two reductions make the column bundleable without changing any result:

- **Positions are implicit.** Every chromosome in the decode is a gapless 1 kb walk from position
  0, so a window's coordinates follow from its row: `[i × 1000, i × 1000 + 1000)`. The header
  records the per-chromosome window counts; the builder asserts the grid property when reading.
- **Values are integers**, `round(P × 10⁶)`. Full float precision costs 5× the space, and the
  probability is never read numerically downstream — it is carried through the beds as a label.

The header (`#` lines) names the individual, the scaling and the per-chromosome counts, so each
file is readable on its own:

```bash
zcat data/YRI_NA18488.PO_5589889.txt.gz | head -30
```

## Building the segments

```bash
python build_human_uncoalesced.py                 # -> ../expected_uncoal/human/
python build_human_uncoalesced.py --skip-phylop   # segments only, no bigWig needed
```

    sum P over the 7 individuals per window  ->  keep sum > 1  ->  drop windows touching a
    segmental duplication  ->  mean phyloP  ->  BMAP  ->  the four "unconstrained" cutoffs


Outputs, written into `../expected_uncoal/human/`, are what `../build_chimp_and_intersect.py`,
`../phylop_annotate_and_gene_test.sh` and the Figure 4 and 6 scripts read:

| file | rows |
|---|---|
| `expected_uncoal_human_5589889.chr{1..22}.bed` | 65,128 total |
| `phylop/…chr{N}.phylop100way.bed`, `…bmap.bed` | + mean phyloP, + BMAP |
| `phylop/…genome_wide.phylop100way.bmap.bed` | 65,128 |
| `phylop/…genome_wide.phylop100way.unconstrained_cutoff_0.331998.bed` | 62,174 |

The other three cutoffs (0.090083, 0.150087, 0.617324) are written too; 0.331998 is the one the
paper reports. Chromosomes are concatenated in string-sorted order (chr1, chr10, chr11, …, chr9), which is the
order the downstream scripts expect.

## Not bundled

`../resources/chm13v2.0_SD.pos_only.hg38_coords_only.bed`, `../resources/hg38.phyloP100way.bw`
and `../resources/whole_genome.1kb.bmap_deciles.bed`. See [`../../../RESOURCES.md`](../../../RESOURCES.md).

## Software

Python 3 (numpy, pandas); bedtools; UCSC `bigWigAverageOverBed`.
