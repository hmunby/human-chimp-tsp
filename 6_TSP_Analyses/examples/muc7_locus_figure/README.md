# MUC7 locus figure (Figure 5)

The MUC7 region is the worked example of a multi-TSP locus: five TSPs across 822 bp
(`4_71210148_71210970`), all on one haplotype (minimum pairwise r² 0.987). Figure 5 shows it
three ways.

| panel | what | made by |
|---|---|---|
| A | one posterior local genealogy at two of the TSPs | `plot_panel_a.py` |
| B | per-SNP TSP support across the locus, on a LocusZoom plot | `build_panel_b_input.py` + the LocusZoom site |
| C | allele frequencies of chr4:71,210,336 by population | not scripted; read from the TSP set |

The three panels are assembled into the final figure by hand.

## Panel A — the local genealogy

The reported panel is **seed 3, sample 134** — a 51 bp tree, TSP branch at 7.7 My, TMRCA
15.5 My. 

```bash
python plot_panel_a.py --seed 3 --sample 134            # -> panelA_seed3_sample134.svg
```

### Allele polarity

SINGER writes states as 0/1, meaning REF/ALT. `--polarity epo` (default) labels mutations
ancestral-to-derived using the EPO ancestral allele, matching how the manuscript reports them;
`--polarity singer` labels the change the tree itself reconstructs along the branch. **The two
are opposite at both sites**, because SINGER places the REF allele at the root while EPO calls
the ALT allele ancestral:

| position | REF | ALT | EPO ancestral | `epo` label | `singer` label |
|---|---|---|---|---|---|
| 71,210,321 | A | T | T | T→A | A→T |
| 71,210,336 | T | C | C | C→T | T→C |

## Panel B — TSP support across the locus

`build_panel_b_input.py` writes the per-SNP pooled support (mean over the three chains, i.e.
the percentage of all 300 posterior genealogies meeting the topology criteria) into the
`neg_log_pvalue` column, so the LocusZoom y-axis reads directly as support (%).

```bash
python build_panel_b_input.py
```

| file | contents |
|---|---|
| `panelB_support.allshared.tsv` | all 35 shared SNPs in the 775 kb SINGER region |
| `panelB_support.allshared.tsv.gz` + `.tbi` | the same, bgzipped and tabix-indexed for LocalZoom |
| `panelB_support.candidates_window.tsv` | the 11 age ≥ 4 Mya candidates within the plotted window |

**Support is not a p-value.** The `pvalue` column is a convenience for tools that demand one.

LocalZoom column mapping (region queries return no header, so specify by number):
chromosome 1, position 2, ref 3, alt 4, p-value 5 with **"-log10 p-value"** ticked. Build
**GRCh37**; chromosome names are unprefixed. Region `4:71,150,000-71,350,000` is the published
window, with `4:71210336` as the LD reference variant.

## Panel C

Frequencies of chr4:71,210,336 in the 1000 Genomes continental groups and the chimpanzee
subspecies. Not scripted: the human and chimpanzee frequencies are the `AF` and `CHIMP_AF`
columns of the TSP set
(`../../../4_SINGER_GenealogicalReconstruction_TSP/4_TSPRegions/2_ClusteringLD/tsps.pooled.hg19.txt`),
and the per-population breakdown is computed from the 1000 Genomes and chimpanzee VCFs.

## Software

Python 3 (tskit ≥ 1.0, pandas, numpy); bcftools; htslib (`bgzip`, `tabix`); `zstd` for the ARG
archives. Rendering SVG to PNG for previews used `cairosvg`; the SVGs themselves need nothing.
