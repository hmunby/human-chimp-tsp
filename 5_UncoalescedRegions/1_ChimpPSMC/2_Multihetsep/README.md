# 1_ChimpPSMC / 2_Multihetsep — alignments to multihetsep

Turns each chimpanzee BAM/CRAM into the per-chromosome multihetsep (`.mhs`) files cobraa reads.

## Rules, in order
1. **get_coverage** — mean depth on chromosome 20, per sample. `bamCaller.py` uses it to set the
   per-site minimum and maximum depth cutoffs, so it must be computed before calling.
2. **convert_bam_to_vcf_and_bed** — `bcftools mpileup -q 20 -Q 20 -C 50` → `bcftools call -c -V indels`
   → `bamCaller.py`, giving a per-chromosome VCF of variant sites **and** the callability BED of
   sites that passed the depth filters. Both are needed: the multihetsep records the number of
   callable bases between consecutive heterozygous sites, so a variant file alone is not enough.
3. **convert_vcfs_beds_to_mhs** (and `_stringent`) — `generate_multihetsep.py` intersects the
   callability BED with the mappability mask and writes the `.mhs`.

Two outputs per sample per chromosome, differing only in mappability stringency:

| output | mask | consumed by |
|---|---|---|
| `{pop}/{sample}/mhs/chr{c}.mhs` | r = 0.50 | `../3_Cobraa/` rule `run_cobraa` (main) |
| `{pop}/{sample}/mhs/stringent75/chr{c}.mhs` | r = 0.75 | `../3_Cobraa/` rule `run_cobraa_high_stringency` |

## Files
- `Snakefile` — the pipeline.
- `config.yaml` — Snakemake cluster (SLURM) profile; edit for your scheduler.
- `misc_files/`
  - `sample_subspecies.txt` — every sample's subspecies (WES/EAS/CEN).
  - `unrelated_samples.txt` — the unrelated individuals to keep.
  - `chr_chimp.txt` — GenBank contig name ↔ chimpanzee chromosome name.

## Inputs to set
Edit the `#### INPUT PATHS ####` block at the top of the `Snakefile`:

- `REFERENCE_FA` — panTro6 / Clint_PTRv2 (GCA_002880755.3).
- `MAPPABILITY_MASK`, `MAPPABILITY_MASK_STRINGENT` — the two BEDs from `../1_MappabilityMask/`.
- `BAM_CALLER`, `GENERATE_MULTIHETSEP` — see below.
- alignments symlinked to `bam/{subspecies}/{sample}.bam|.cram` with indexes. Samples stored as
  CRAM are listed in `SAMPLES_CRAM`.

### msmc-tools (not bundled)
`bamCaller.py` (with its `utils.py`) and `generate_multihetsep.py` come from
[msmc-tools](https://github.com/stschiff/msmc-tools). Place them in `scripts/`.

`generate_multihetsep.py` needs one local change: **an added `--chr` argument** that overrides the
chromosome name written to the output. The panTro6 reference uses GenBank contig names
(`CM009238.2`), and cobraa and everything downstream expect chimpanzee chromosome names (`1`,
`2A`, …), so the rules pass `--chr={wildcards.chrom}` to rename on the fly. Without that patch the
`.mhs` files carry contig accessions and the decode output cannot be lifted.

## Software
bcftools 1.20; samtools; msmc-tools; Python 3; Snakemake 7.
