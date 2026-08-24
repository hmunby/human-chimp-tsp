# 5 — Uncoalesced regions

Genomic regions that have not coalesced between human and chimpanzee by 5.5 Mya, from the
chimpanzee cobraa/PSMC decode, and their reciprocal (human × chimp) overlap, phyloP
conservation, BMAP dependence and gene-set enrichment.

## Pipeline

### 1_ChimpPSMC — alignments → multihetsep → cobraa → decode → minimal decode → hg38

The chimpanzee inference, in four sub-stages (see its own README; "PSMC" here means **cobraa**,
which takes multihetsep rather than `.psmcfa` input):

- `1_MappabilityMask/` — SNPable k=150 panTro6 mask at r=0.50 (main) and r=0.75 (sensitivity).
- `2_Multihetsep/` — BAM/CRAM → `bcftools mpileup`/`call` → `bamCaller.py` → per-chrom VCF +
  callability BED → `generate_multihetsep.py` → `.mhs`.
- `3_Cobraa/` — cobraa fit (D=64, b=100, spread1=0.075, spread2=100, muoverr=1.5, 30 iters,
  thetafixed=0.001) for 9 chimpanzees (CEN/EAS/WES), then the posterior decode → per-sample,
  per-chrom posterior grid (panTro6).
- `4_MinimalDecode/` — mu=1.27e-8 → N=19685 → grid index 49 → age 5,755,625 yr → selects the
  P(uncoalesced beyond 5.5 Mya) column (`PO_5755625`) → minimal-decode beds (panTro6); then
  liftOver panTro6→hg38, keeping 950–1050 bp segments → `cleaned_lift_hg38/`.

### 2_UncoalescedSegments — both species' segments, reciprocal intersect, phyloP and genes

- `human_decode/` — the bundled Cousins et al. decode column for the 7 humans, and
  `build_human_uncoalesced.py`: sum P(uncoal) across the 7 per 1 kb (keep sum > 1), drop
  segmental duplications, annotate phyloP and BMAP → human genome-wide uncoalesced bed and the
  "unconstrained" cutoffs.
- `build_chimp_and_intersect.py` — sum P(uncoal) across the 9 samples per 1 kb (keep sum > 1),
  merge, annotate BMAP → chimpanzee genome-wide uncoalesced bed; reciprocal-intersect with the
  human uncoalesced bed.
- `merge_bed_weighted.py` — length-weighted bed merge, used to build the closest-gene test set.
  Merges only strictly bookended segments, which `bedtools merge` does not.
- `phylop_annotate_and_gene_test.sh` — two independent branches. **(A)**
  `bigWigAverageOverBed` → phyloP < 0.331998 → reciprocal unconstrained intersect → **Figure 6**.
  **(B)** length-weighted merge of the *unfiltered* exact intersect → `bedtools closest`
  (GENCODE v49) → the enrichment input.

### 3_Figures

- `reciprocal_uncoal_deciles.py` — reciprocal uncoalescence by BMAP decile (**Figure 4A**), and
  the TSP counts and assigned genes per decile (**4B**, **4C**).
- `plot_fig6_obs_vs_exp.py` — observed versus expected reciprocal uncoalescence by BMAP quartile
  (**Figure 6**); `.withMHC.py` is the MHC-included supplementary version.
- `per_individual/` — percentage of each genome uncoalesced by 5.5 Mya, per individual and
  decile (**Supplementary Figure 5**).

The overlap of the SNP sets (TSP, old shared, all shared, 1000GP) with these regions is
**Figure 2**, and lives with the other SNP-set comparisons at
`../6_TSP_Analyses/snp_set_comparisons/uncoalesced_overlap.py`. It reads the reciprocal
intersect produced here.

### 4_Enrichment

Gene enrichment among the genes nearest the uncoalesced regions, per BMAP quartile. Two tests,
as reported:

- `enrichment_uncoal_naive_and_perm.py` — the **naive hypergeometric**: test genes against all
  genes nearest any callable segment in the same quartile. This is the test that shows
  glycoproteins enriched in every quartile.
- `enrichment_uncoal_locus_sweep.py` — the **matched permutation**, at locus level. Regions are
  merged into loci at distance D, one gene per locus, and compared to null loci drawn to match
  span within the same quartile, so that long and isolated genes are drawn as often as they are
  captured. Swept over D and over a distance cap.
- `glyco_genelevel_hypergeom.py` and `make_glyco_summary_table.py` — the glycoprotein and
  membrane-glycoprotein tables, joining the two tests for the same conditions.
- `collapse_redundant_categories.py` — gene-set libraries are highly redundant, so significant
  categories are grouped by the genes actually driving them (EnrichmentMap combined coefficient
  0.5·Jaccard + 0.5·overlap ≥ 0.375).
- `make_full_results_tables.py` — one full-results table per merge distance: every significant
  category, its driver genes and its cluster.

## Not bundled

The BMAP decile bed, the hg38 phyloP bigWig, and the uncoalesced-segment and all-segment BEDs. Place them under each sub-stage's `resources/`.

## Software

cobraa; bedtools v2.29.1; UCSC `liftOver` and `bigWigAverageOverBed`; Python 3 (numpy, pandas,
matplotlib, scipy); Snakemake 7.
