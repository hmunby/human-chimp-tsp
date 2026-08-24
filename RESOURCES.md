# External inputs (`resources/`)

Every script in this repository reads its large inputs from a `resources/` directory beside it, so
that nothing points outside the repository. Those inputs are not bundled: they are public
reference data, or intermediates that earlier stages produce, and together they run to several
hundred gigabytes.

Each stage directory that needs them expects its own `resources/`. The simplest arrangement is one
real directory somewhere with space, symlinked into each stage:

```bash
mkdir -p /path/to/resources
ln -s /path/to/resources 5_UncoalescedRegions/2_UncoalescedSegments/resources
```

`common/config.py` takes `PAPER_DATA_ROOT` as an override if you would rather keep one root and
point at it:

```bash
export PAPER_DATA_ROOT=/path/to/resources
```

Anything the scripts *write* also lands under `resources/` , so it needs to be writable, not read-only.

---

## 1. Public reference data

Download once; shared by several stages.

| what | where it goes | source |
|---|---|---|
| panTro6 / Clint_PTRv2 FASTA | `resources/GCA_002880755.3_Clint_PTRv2_genomic.fa`, `resources/panTro6.simple_chrom.no_chr.fa` | GenBank GCA_002880755.3 |
| hg19 FASTA + chrom sizes | `resources/hg19/ucsc.hg19.fasta`, `resources/hg19/human.hg19.genome` | UCSC |
| hg38 FASTA + chrom sizes | `resources/hg38/hg38.fa`, `resources/hg38.genome` | UCSC |
| liftOver chains | `resources/chains/` — `panTro6ToHg19.no_chr.over.chain`, `hg19ToPanTro6.over.chain`, `hg19ToHg38.over.chain`, `b37tohg38.nochr.over.chain`, `panTro6ToHg38.over.chain` | UCSC |
| 1000 Genomes phase 3 VCFs + sample panel | `resources/1000GP/` | EBI, 20130502 release |
| GENCODE v49 annotation | `resources/gencode.v49.annotation.gtf.gz` and the derived beds below | GENCODE |
| hg38 phyloP 100-way bigWig | `resources/hg38.phyloP100way.bw` | UCSC |
| phyloP 36-mammal EPO, human excluded | `resources/phylop_36mammal_no_human/…chr{1..22}.phyloP_no_human.txt` | Fu et al. 2014 |
| phyloP 46-way placental | `resources/phylop_46way/placentalMammals.phyloP46way.chr{1..22}.txt` | UCSC |
| ENCODE SCREEN cCREs | `resources/functional_annotation_datasets/` | SCREEN V4 |
| B-map (CADD best-fit) | `resources/whole_genome.bmap.bed.gz` | Murphy et al. 2023 |
| Roulette mutation rates | `resources/roulette/autosomes/` | Seplyarskiy et al. 2023 |
| HapMap II + interpolated OMNI genetic maps | `resources/1000-genomes-genetic-maps/{hapmapII,interpolated_OMNI}/` | github.com/joepickrell/1000-genomes-genetic-maps |
| T2T segmental duplications | `resources/chm13v2.0_SD.pos_only.hg38_coords_only.bed` | Vollger et al. 2022 |
| snpEff 5.1 + hg19 database | `resources/snpEff/snpEff.jar` | snpEff |

**SINGER posterior allele ages.** Deng, Nielsen & Song (2025, *Nat Genet* 57:2124–2135), 100
posterior ARG samples over 200 African 1000GP haplotypes, used by the supplementary age figures.
On Zenodo as per-chromosome tarballs: [10437053](https://zenodo.org/records/10437053) (chr1–4),
[10467284](https://zenodo.org/records/10467284) (chr5–10, 22),
[10467509](https://zenodo.org/records/10467509) (chr11–21). Unpack and reduce to the age pickles at
`resources/ages_all/african_chr{N}_ages.pickle`; they need numpy 2.x to unpickle. See
`3_SharedVariation/figures/supplementary/README.md`.

**Shared-SNP annotation table** (`all_shared_snps_ann.txt`). One row per shared SNP in hg19,
carrying REF/ALT, the Ensembl ancestral allele, human and chimpanzee allele frequencies, BMAP,
CpG status, the Roulette mutation rate and the Relate age bounds. Produced by stage 3 and read by
the TSP calling, Figure 1 panel B and the phyloP confound test.

**Derived GENCODE beds.** Several scripts want protein-coding autosomal genes as BED with gene
names. Build once from the GTF and place all of:
`resources/gencode.v49.annotation.genes.protein_coding.autosomes.gene_names.bed`,
`resources/gencode.v49.pc_genes.bed`, `resources/gencode.v49.genes.no_readthrough.bed`,
`resources/gencode.v49.genes.protein_coding.autosomes.no_readthrough.bed`.

## 2. Software on `PATH`

`bwa`, `samtools`, `bcftools` 1.20 (with the `+liftover` plugin), `bedtools` 2.29.1,
`PLINK` 1.90, UCSC `liftOver` / `bigWigAverageOverBed` / `bigWigInfo`, `snakemake` 7, `R` with
tidyverse, and Python 3 with numpy / pandas / scipy / matplotlib / biopython / msprime / tskit.

Three tools are cloned rather than installed, and their paths are set at the top of the scripts
that use them: **cobraa** (`resources/cobraa/cobraa.py`), **SINGER**, and **msmc-tools**
(`bamCaller.py`, `utils.py`, `generate_multihetsep.py` into
`5_UncoalescedRegions/1_ChimpPSMC/2_Multihetsep/scripts/`; note the local `--chr` patch described
in that directory's README). **SNPable** is needed for the mappability mask.

## 3. Stage outputs consumed by later stages

These are produced by running the pipeline in order. If you are re-running only part of it, these
are what the later stage needs from the earlier one.

| produced by | needed as |
|---|---|
| stage 1 filtering | `resources/genotyped.chimp43.vcf.gz` (panTro6), `resources/genotyped.chimp43.hg19.b37.vcf.gz` |
| stage 2 annotation | `resources/1000GP.wgs.biallelic.sites.YRI.CpG.snpeff.MR.hg38.MAF05.vcf.gz`, `resources/1000GP.wgs.biallelic.sites.MAF05.vcf.gz` |
| stage 3 intersect | `resources/shared_variants.SD_CNV_DUP_filtered.vcf.gz`, `…sites.txt`, `…snpeff.txt`, `resources/shared_snps_chrom_pos_hg38.sorted.bed`, `resources/all_shared_snps_ann.txt` |
| stage 4 chimp recall | `resources/chimp_bam_hg19/`, phased haplotypes |
| stage 5 `1_ChimpPSMC/4_MinimalDecode` | `resources/cleaned_lift_hg38/` — the per-sample hg38 minimal-decode beds |
| stage 5 `2_UncoalescedSegments` | `resources/expected_uncoal/`, both species — `chimp/` and the reciprocal intersect from `build_chimp_and_intersect.py`, `human/` from `human_decode/build_human_uncoalesced.py`. Both halves must exist before `phylop_annotate_and_gene_test.sh` will run. Read by the stage 5 figures and enrichment, and by `6_TSP_Analyses/snp_set_comparisons/uncoalesced_overlap.py` |
| stage 5 mappability mask | `resources/PanTro6_mappability_mask.150.{50,75}.bed` |
