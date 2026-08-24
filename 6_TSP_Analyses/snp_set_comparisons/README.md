# snp_set_comparisons — the TSPs against the other SNP sets

Every analysis here contrasts the final TSP set with the progressively larger SNP sets it is
drawn from, so the same four nested sets and the same palette recur throughout:

| set | n | colour | definition |
|---|---|---|---|
| 1000GP (MAF ≥ 5%) | 6,435,348 | grey `#7F7F7F` | 1000 Genomes biallelic autosomal sites, MAF ≥ 5% |
| All shared | 66,622 | blue `#1F77B4` | human–chimp shared SNPs (stage 3) |
| Old shared (≥ 4 Mya) | 860 | purple `#9467BD` | shared SNPs with YRI midpoint age ≥ 4 Mya, MHC excluded (stage 3 candidates) |
| TSP | 117 | red `#FF0000` | the reported (pooled) trans-species set (stage 4) |

A fifth track, **chimpanzee (AC ≥ 3)**, appears in the phyloP figures only (black); it has no
counterpart in the variant-effects panel.

Because old-shared and TSP are *subsets* of the shared SNPs, they are never scored separately:
their positions are looked up in the shared table. That guarantees identical treatment and is
why the pipeline below only annotates three sets.

## Components

### Building the phyloP inputs
- `Snakefile` — annotates the shared / 1000GP-MAF5 / chimpanzee SNP positions (hg19) with
  phyloP, against two tracks: `phylop/` (36-eutherian-mammal EPO, **human sequence excluded** —
  the main-figure track) and `phylop_46way/` (UCSC 46-way placental, human included). Output is
  one table per set per chromosome, which every phyloP script below reads.
- `scripts/annotate_phylop_scores.py` — the per-chromosome left join.
- `config.yaml` — Snakemake cluster (SLURM) profile; edit for your scheduler.

### Analyses
- `uncoalesced_overlap.py` — **Figure 2**: percentage of each set falling in the reciprocal
  (human × chimp) uncoalesced segments from stage 5, against a size-matched
  1000GP null (1000 permutations).

  | set | n | in uncoalesced | fold | p |
  |---|---|---|---|---|
  | 1000GP (MAF ≥ 5%) | 6,435,348 | 0.736% | — (null pool) | — |
  | All shared | 66,622 | 1.53% | 2.08× | 9.99e-4 |
  | Old shared (≥ 4 Mya) | 858 | 16.6% | 22.9× | 9.99e-4 |
  | TSP | 117 | 47.0% | 66.7× | 9.99e-4 |

  Outputs are bundled: `snp_set_uncoal_overlap.{txt,png,no_expected.png}`.
- `variant_effects.py` — **Panel A**: snpEff consequence distribution per set, with the
  greatest-consequence rule re-derived per SNP from the raw multi-annotation field (not
  snpEff's four-level IMPACT). Categories ordered by the Agarwal et al. Fig 2B constraint
  estimate. Produces full-category and merged-category versions.
- `ccre_overlap.py` — percentage of each set overlapping an ENCODE SCREEN candidate
  cis-regulatory element (GRCh38). The reported registry is **V4** (Moore et al., *Nature*
  2026; 2,348,854 elements, 21.2 % of the autosomes); `--registry v3` scores against the 2020
  registry instead.

  | set | overlapping a cCRE |
  |---|---|
  | 1000GP (MAF ≥ 5%) | 22.3% |
  | All shared | 20.7% |
  | Old shared (≥ 4 Mya) | 18.5% |
  | TSP (pooled, 117) | 15.4% |

  Output: `snp_set_ccre_overlap.txt`.
- `phylop_ecdf.py` — **Panel B**: five-track phyloP eCDF, for both phyloP tracks.
- `combined_variant_effects_phylop.py` — assembles Panels A and B into the combined figure
  (side-by-side and stacked), importing the two scripts above so the panels cannot drift from
  the standalone versions.
- `phylop_stats.py` — the descriptive statistics and tests behind the phyloP figure: two-sided
  Mann-Whitney with effect sizes plus a KS test per pair, BH-corrected, reported with and
  without the MHC and for both tracks. Redirect stdout to `phylop_stats_output.txt`.
- `phylop_confound_test.py` — supplementary control: separates the chimpanzee
  reference-allele artefact from CpG hypermutability as explanations for the shared-SNP phyloP
  deficit, by stratifying on chimp derived-allele frequency and CpG status.

## Run order

```bash
snakemake --profile .            # phylop/ + phylop_46way/ score tables
python variant_effects.py        # Panel A
python phylop_ecdf.py            # Panel B
python combined_variant_effects_phylop.py
python phylop_stats.py > phylop_stats_output.txt
python uncoalesced_overlap.py    # Figure 2
python ccre_overlap.py           # SCREEN cCRE overlap (V4)
python phylop_confound_test.py   # supplementary
```

`variant_effects.py` and `phylop_ecdf.py` build their data at import, so the combined-figure
driver needs no separate configuration.

## Bundled inputs (`misc_files/`)

| file | rows | used by |
|---|---|---|
| `tsp_snps_pooled.hg19.txt` | 117 | variant effects, phyloP eCDF/stats/confound |
| `tsp_snps_pooled.hg38.bed` | 117 | uncoalesced overlap |
| `tsp_snps_{strict,ne20k,ne100k}.{hg19.txt,hg38.bed}` | 85 / 111 / 133 | the sensitivity sets, selectable with `--tsp-file` / `--tsp-bed` |
| `old_shared_snps_860.hg19.chrom_pos.txt` | 860 | variant effects, phyloP eCDF/stats |
| `old_shared_snps_860.hg38.sorted.bed` | 858 | uncoalesced overlap (the 860 set lifted; see below) |
| `agarwal_et_al_cpg_fig2b.txt` | 7 | consequence ordering (constraint estimate per category) |
| `tgp_annotation_simple_counts.txt` | 12 | precomputed 1000GP category counts, built with the same greatest-consequence rule over its 2.8 GB snpEff file |

## Inputs to set (`resources/`)

Not bundled — point these at your local copies, or at the archives named in the data-availability
statement:

- the annotated shared VCF and the 1000GP MAF ≥ 5% VCF (stages 3 and 2), plus the hg19-lifted
  chimpanzee VCF (stage 1) — the `Snakefile`'s three SNP sources;
- the two phyloP tracks, one file per chromosome (`CHROM<TAB>POS<TAB>SCORE`);
- `shared_variants.SD_CNV_DUP_filtered.snpeff.txt` (~27 MB) — `bcftools query -f '%CHROM\t%POS\t%ANN\t%LOF\n'`
  over the annotated shared VCF, for `variant_effects.py`;
- `tgp_snps_hg38.bed` (153 MB) and `shared_snps_chrom_pos_hg38.sorted.bed` (1.6 MB) — the two
  large SNP beds for `uncoalesced_overlap.py`;
- `expected_uncoal_intersect.whole_genome.exact.bed` — the reciprocal intersect
  from `../../5_UncoalescedRegions/2_UncoalescedSegments/`;
- `all_shared_snps_ann.txt` — the shared-SNP annotation table (hg19), for
  `phylop_confound_test.py`.

## The old-shared set is 860 everywhere

"Old shared" means the **860** stage-3 candidate SNPs: shared SNPs with a YRI midpoint Relate age
of at least 4 Mya, MHC excluded. Every analysis in this directory uses that set.

In hg38 the set is **858** SNPs. Two do not lift and are dropped: `9:136143120` and
`14:106529056`. This is liftOver attrition only, not a different definition.

The candidate definition is applied once, in `../../3_SharedVariation/2_CandidateTSPs/`; nothing
downstream re-filters it. Of the 858, **142** fall in the reciprocal uncoalesced segments (16.6%,
22.9× the size-matched 1000GP expectation).

## Software

Python 3 (numpy, pandas, matplotlib, scipy); bedtools v2.29.1; bcftools 1.20; Snakemake 7.
