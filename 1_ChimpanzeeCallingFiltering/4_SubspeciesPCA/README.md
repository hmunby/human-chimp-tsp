# 1 / 4_SubspeciesPCA — chimpanzee subspecies identification (Supplementary Figure 1)

The publicly available chimpanzees carry subspecies labels from their source publications; the five
founders of the newly sequenced pedigree do not. This assigns them by clustering.

PCA on LD-pruned genome-wide genotypes separates the unrelated samples into three clusters
corresponding to the three subspecies. **PC1 explains 33.8% of the variance and PC2 11.7%.** Four of
the five new founders fall in the Western cluster and one in the Eastern cluster, giving the final
set of 14 Eastern, 14 Central and 15 Western chimpanzees.

## Flow (rules, in order)
1. **rename_chroms** — contig accessions to chromosome names.
2. **plink_make_bed** — VCF to plink binary.
3. **plink_linkage_prune** — `--indep-pairwise 50 10 0.1`.
4. **plink_pca** — PCA on the pruned set.
5. **get_subspecies_ids** — assign each sample to a subspecies from its PC1/PC2 position.
6. **plot_pca** — the figure.

Cluster boundaries, read off the already-labelled samples: WES at PC1 > 0; CEN at PC1 < 0 and
PC2 < 0; EAS at PC1 < 0 and PC2 > 0.

## Files
- `Snakefile` — the pipeline.
- `get_subspecies.R` — cluster-to-subspecies assignment.
- `plot_pca.R` — Supplementary Figure 1: PC1 vs PC2, coloured by subspecies, shaped by whether the
  label was already known (circle) or assigned here (triangle).
- `misc_files/`
  - `rename_chroms.txt`, `sample_known_subspecies.txt` (`UNK` marks the unlabelled founders).
  - `pca.eigenvec`, `pca.eigenval` — **the PCA output as run**, bundled so the figure can be
    regenerated without the 4.6 GB VCF:
    ```
    Rscript plot_pca.R misc_files/pca.eigenvec misc_files/pca.eigenval \
                       misc_files/sample_known_subspecies.txt pca.pdf
    ```

## Inputs to set
`VCF` at the top of the `Snakefile` — the combined filtered chimpanzee VCF from `../3_Filtering/`.

## Software
PLINK v1.90; bcftools 1.20; R with tidyverse; Snakemake 7.
