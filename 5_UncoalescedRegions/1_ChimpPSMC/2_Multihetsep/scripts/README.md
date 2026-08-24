# scripts/ — msmc-tools

Place `bamCaller.py`, its `utils.py`, and `generate_multihetsep.py` here, from
[msmc-tools](https://github.com/stschiff/msmc-tools). They are not bundled.

`generate_multihetsep.py` needs one local change: **add a `--chr` argument** that overrides the
chromosome name written to the output. The panTro6 reference uses GenBank contig names
(`CM009238.2`), while cobraa and everything downstream expect chimpanzee chromosome names
(`1`, `2A`, …), so the Snakefile passes `--chr={wildcards.chrom}` to rename on the fly.

Without that patch the `.mhs` files carry contig accessions and the decode output cannot be
lifted to hg38.
