# misc_files — bundled inputs for the intersection

- `AF_ID_header.hdr` — INFO header lines for the chimp annotations added to the shared-variant VCF
  (`CHIMP_AF` = chimp allele frequency; `CHIMP_ID` = SNP position in the panTro6 assembly).
- `paralogs_SD.bed` — **human** segmental duplications: the CHM13v2.0 (T2T) SD annotation, lifted to
  hg19 coordinates.
- `paralogs_CNV_DUP.bed` — **human** duplication / CNV regions: DUP + CNV records from the 1000
  Genomes phase 3 structural-variant callset (`ALL.wgs.mergedSV.v8.20130502`).

These are **human** paralog masks (species A side), applied to remove paralog-driven false shared
SNPs. They are distinct from the chimpanzee segmental-duplication track used in stage 1.

Merged / union versions (used only for total-length counting, not by the filter rules) can be
regenerated with bedtools:
```
bedtools merge -i paralogs_SD.bed > paralogs_SD.merged.bed
bedtools sort -i paralogs_CNV_DUP.bed | bedtools merge -i - > paralogs_CNV_DUP.merged.bed
cat paralogs_SD.merged.bed paralogs_CNV_DUP.merged.bed | bedtools sort -i - | bedtools merge -i - > paralogs_SD_CNV_DUP.merged.bed
```
