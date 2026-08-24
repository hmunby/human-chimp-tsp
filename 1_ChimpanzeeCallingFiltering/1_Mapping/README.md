# 1_Mapping

Align chimpanzee reads to panTro6 (GCA_002880755.3) and mark duplicates.

1. `1_bwa_mem.sh` — `bwa mem` (read groups added inline via `-R`) piped to `samtools sort`, one
   task per sequencing run. Run once per project (set `PROJECT` + its sample/run ID lists).
2. (per-sample merge) — where a sample has more than one run, merge its per-run BAMs
   (`samtools merge`) into one BAM before duplicate marking.
3. `2_mark_duplicates.sh` — GATK/Picard `MarkDuplicates` per sample → deduplicated BAMs (the input
   to `../2_VariantCalling/`).


Paths and sample lists are set in the parameter block at the top of each script.
