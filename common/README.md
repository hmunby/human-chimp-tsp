# common — shared figure infrastructure

Shared modules imported by the figure scripts across the repo.

- `config.py` — central registry of input-data paths + shared constants (generation time, the
  human–chimp split band, Roulette scaling, etc.). Paths default to the analysis tree and are
  overridable without editing:
  ```bash
  export PAPER_DATA_ROOT=/path/to/chimp_human_shared_variation
  export PAPER_FIG_OUTDIR=/path/to/output
  ```
- `plotting.py` — shared plotting helpers

Figure scripts live with their stage/analysis and import these; run them with `common/` on
`PYTHONPATH` (e.g. `export PYTHONPATH=$PWD/common`), or from a directory where `import config` /
`import plotting` resolve. Requires `numpy`, `pandas`, `matplotlib`.
