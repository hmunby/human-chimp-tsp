# 3 — Shared human–chimpanzee variation

Builds the trans-species shared-SNP set from the human (stage 2) and chimpanzee (stage 1) call sets,
and the analyses that follow from it.

## Components
- `1_Intersect/` — intersect human × chimp SNPs, filter human paralogs, and annotate → the
  shared-variant set.
- `2_CandidateTSPs/` — identify candidate TSPs (>=4 Mya in Relate YRI)
