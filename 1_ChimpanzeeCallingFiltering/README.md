# 1 — Chimpanzee variant calling and filtering

Stage 1 of the pipeline: from public chimpanzee FASTQ to a filtered, hg19-lifted chimp variantset (the chimp side of the human–chimp shared-SNP intersection used throughout the paper).

Reference: **panTro6** (GenBank GCA_002880755.3, "Clint_PTRv2"). Public chimpanzee short-read data come from several projects: PRJEB5937, PRJEB29710, PRJDB3537, PRJEB15086, and PRJNA1357362 (10 individuals first published in this study). The panel is 59 individuals, 43 remain after removing related individuals (see the filtering step).

## Pipeline
```
FASTQ ─(1_Mapping)─▶ deduplicated BAMs ─(2_VariantCalling)─▶ raw joint-genotyped VCFs
      ─(3_Filtering, Snakemake)─▶ filtered VCF ─▶ lifted panTro6→hg19  (input to the intersection)
```

- **1_Mapping/** — `bwa mem` (read groups added inline) → `samtools sort` → per-sample merge →
  GATK/Picard `MarkDuplicates`.
- **2_VariantCalling/** — GATK `HaplotypeCaller` (GVCF) → `GenomicsDBImport` → `GenotypeGVCFs`.
- **3_Filtering/** — the VCF-processing Snakemake: biallelic → site/genotype QC →
  Mendelian-consistency filter → subset to unrelated individuals → segmental-duplication filter →
  allele-count filter → combine → liftOver to hg19.
- **4_SubspeciesPCA** - PCA to determine subspecies of individuals published in this study.  

## Supporting inputs 
FASTQ were obtained from the listed ENA/DDBJ projects and the reference from GenBank (GCA_002880755.3), liftOver chains (panTro6↔hg19) from UCSC. 
