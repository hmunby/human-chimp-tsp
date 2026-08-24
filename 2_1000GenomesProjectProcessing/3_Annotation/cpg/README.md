# CpG annotation

Builds `cpg_annot_table.txt.gz` (`CHROM,POS,REF,ALT,CPG`), used by the annotation driver to add the
`CPG` INFO field. `CPG = 1` marks a C→T transition in a CpG dinucleotide context.

CpG status is determined from the **ancestral** allele (1000GP `AA` field)

## Steps
1. `make_cpg_annotation_table.sh` (SLURM array 1–22) → per-chrom `cpg_tables_per_chrom/cpg.annot_table.chr{N}.txt`.
   Each chromosome is scored by `make_cpg_annotation.R`.
2. `combine_cpg_tables.sh` → concatenate chr1–22, `bgzip`, `tabix -s1 -b2 -e2` → `cpg_annot_table.txt.gz`.

## Software
- R 4.5.1, Bioconductor 3.21: `MutationalPatterns`,
  `VariantAnnotation`, `BSgenome`, `BSgenome.Hsapiens.1000genomes.hs37d5`, `stringr`.
- bcftools 1.20 (htslib 1.20); tabix.

The CpG dinucleotide context is read from the `BSgenome.Hsapiens.1000genomes.hs37d5` (GRCh37/hs37d5)
reference; the ancestral allele is the 1000GP `AA` INFO field (Ensembl EPO-derived).
