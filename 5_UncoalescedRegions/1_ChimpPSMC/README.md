# 5 / 1_ChimpPSMC — chimpanzee demographic inference and posterior decoding

From chimpanzee alignments to the per-base probability that an individual's two haplotypes have
not yet coalesced by the human–chimp split. Those probabilities are what
`../2_UncoalescedSegments/` turns into the chimpanzee uncoalesced regions.

## Flow

```
BAM/CRAM ──1──> mappability mask (SNPable)
         └──2──> per-chrom VCF + callability BED ──> multihetsep (.mhs)
                        └──3──> cobraa fit (per individual) ──> posterior decode
                                       └──4──> P(uncoalesced at the split) ──> liftOver hg38
```

1. **`1_MappabilityMask/`** — SNPable k=150 mask of panTro6 at two stringencies (r = 0.50 for the
   main analysis, r = 0.75 for the sensitivity fit), converted to BED.
2. **`2_Multihetsep/`** — per-sample consensus calling (`bcftools mpileup`/`call` → `bamCaller.py`)
   giving a VCF plus a matched callability BED, then `generate_multihetsep.py` to combine those
   with the mappability mask into `.mhs`.
3. **`3_Cobraa/`** — run PSMC  implemented in cobraa, fit per individual (D=64, b=100, spread₁=0.075, spread₂=100,
   μ/ρ=1.5, 30 iterations, θ fixed at 0.001), then the posterior decode along each chromosome.
4. **`4_MinimalDecode/`** — collapse the posterior grid to the single "not coalesced by 5.5 Mya"
   column, liftOver panTro6 → hg38, and drop segments distorted by the lift.

## Samples

| subspecies | samples |
|---|---|
| CEN (*P. t. troglodytes*) | SAMEA4374772, SAMEA4374778, SAMEA4374785 |
| EAS (*P. t. schweinfurthii*) | 2003, SAMEA4374790, SAMEA4374797 |
| WES (*P. t. verus*) | 1059, SAMEA5204228, SAMEA2421542 |

Reference is panTro6 / Clint_PTRv2 (GCA_002880755.3), whose contigs carry GenBank names — `2_Multihetsep/misc_files/chr_chimp.txt`
maps them to chromosome names.

## Alignments

The BAMs/CRAMs come from `../../1_ChimpanzeeCallingFiltering/1_Mapping/` (bwa mem → MarkDuplicates
against panTro6). `2_Multihetsep/` expects them symlinked into `bam/{subspecies}/{sample}.bam`
(or `.cram`) with indexes alongside.

## Software

cobraa; msmc-tools (`bamCaller.py`, `generate_multihetsep.py`); SNPable (`splitfa`,
`gen_raw_mask.pl`, `gen_mask`); bwa; samtools; bcftools 1.20; UCSC `liftOver`; Python 3 (numpy,
pandas, matplotlib, biopython); Snakemake 7.
