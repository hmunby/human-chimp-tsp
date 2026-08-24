"""
Permutation test for glycoprotein ascertainment bias.

Tests whether the apparent enrichment of glycoprotein genes (KW-0325) in the
uncoalescent region gene set is due to the genomic distribution of the input
regions rather than biology.

Approach:
  1. Merge contiguous input segments.
  2. Compute observed glycoprotein fraction from the real closest-gene output.
  3. Run two sets of permutations:
       a. Unrestricted: shuffle merged segments randomly across autosomes.
       b. SD-excluded:  shuffle but disallow placement in segmentally duplicated regions
          (matching the exclusion applied to the real input).
     Each permutation: bedtools shuffle → bedtools closest → unique genes → glyco fraction.
  4. Empirical p-values = fraction of permutations >= observed.
  5. Save results TSV and comparison plot.

Usage:
  python glycoprotein_permutation_test.py [--n-perms N] [--seed S]
"""

import argparse
import os
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Paths ────────────────────────────────────────────────────────────────────

INPUT_BED = "resources/expected_uncoal/intersect/expected_uncoal_intersect.whole_genome.exact.bed"
CLOSEST_BED = "resources/expected_uncoal/intersect/expected_uncoal_intersect.whole_genome.exact.closest_genes.bed"
GENCODE_BED = "resources/gencode.v49.annotation.genes.protein_coding.autosomes.gene_names.bed"
GENOME_FILE = "resources/hg38.genome"
SD_BED = "resources/chm13v2.0_SD.pos_only.hg38_coords_only.bed"
GLYCO_TSV = "../../../6_TSP_Analyses/enrichment/glycoprotein/human_glycoproteins_KW0325_with_coords.tsv"
OUT_DIR = Path(
    "resources/glycoprotein_distribution"
)

# ── Helpers ──────────────────────────────────────────────────────────────────

def run(cmd: str, **kwargs) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, shell=True, check=True,
                          capture_output=True, text=True, **kwargs)


def glyco_fraction(closest_bed_path: str, glyco_names: set) -> tuple[float, int, int]:
    """Return (fraction, n_glyco, n_total_unique) from a closest-gene BED.

    Gene name is always the second-to-last column (last is the -d distance).
    """
    df = pd.read_csv(closest_bed_path, sep="\t", header=None, dtype=str)
    gene_col     = df.iloc[:, -2]
    unique_genes = set(gene_col.dropna()) - {".", "-1"}
    n_total  = len(unique_genes)
    n_glyco  = len(unique_genes & glyco_names)
    frac     = n_glyco / n_total if n_total > 0 else 0.0
    return frac, n_glyco, n_total


def make_autosome_genome(genome_file: str, tmp_dir: str) -> str:
    """Write a genome file containing only autosomal chromosomes (chr1–chr22)."""
    out = os.path.join(tmp_dir, "autosomes.genome")
    with open(genome_file) as f_in, open(out, "w") as f_out:
        for line in f_in:
            chrom = line.split("\t")[0]
            if chrom.startswith("chr") and chrom[3:].isdigit():
                f_out.write(line)
    return out


def prepare_sd_excl(sd_bed: str, tmp_dir: str) -> str:
    """Add chr-prefix, filter to autosomes, sort and merge SD regions.

    The source file uses bare chromosome numbers (1, 2, … 22) with no chr prefix.
    bedtools shuffle expects names matching the genome file (chr1 … chr22).
    """
    out = os.path.join(tmp_dir, "sd_excl.bed")
    run(
        f"awk '$1+0>=1 && $1+0<=22 {{print \"chr\"$1\"\\t\"$2\"\\t\"$3}}' {sd_bed}"
        f" | sort -k1,1V -k2,2n | bedtools merge -i stdin > {out}"
    )
    result = run(f"wc -l < {out}")
    print(f"SD mask (regions blocked from shuffle placement, merged autosomes): "
          f"{result.stdout.strip()}")
    return out


def merge_bed(bed_in: str, bed_out: str) -> int:
    """Sort, merge, and return number of merged segments."""
    run(f"sort -k1,1V -k2,2n {bed_in} | bedtools merge -i stdin > {bed_out}")
    result = run(f"wc -l < {bed_out}")
    return int(result.stdout.strip())


def closest_genes(a_bed: str, b_bed: str, out_bed: str) -> None:
    run(f"sort -k1,1V -k2,2n {a_bed} | bedtools closest -a stdin -b {b_bed} -d > {out_bed}")


def run_permutations(
    merged_bed: str,
    autosome_genome: str,
    glyco_names: set,
    n_perms: int,
    rng: np.random.Generator,
    label: str,
    excl_bed: str | None = None,
) -> np.ndarray:
    """Run n_perms shuffle → closest cycles and return array of glyco fractions."""
    tmp_shuffled = merged_bed + ".shuffled.bed"
    tmp_closest  = merged_bed + ".closest.bed"

    excl_flag = f"-excl {excl_bed}" if excl_bed else ""
    fracs: list[float] = []

    for i in range(n_perms):
        iter_seed = int(rng.integers(1, 2**31))
        run(
            f"bedtools shuffle -i {merged_bed} -g {autosome_genome} "
            f"{excl_flag} -seed {iter_seed} > {tmp_shuffled}"
        )
        closest_genes(tmp_shuffled, GENCODE_BED, tmp_closest)
        frac, _, _ = glyco_fraction(tmp_closest, glyco_names)
        fracs.append(frac)

        if (i + 1) % 100 == 0:
            print(f"  [{label}] {i+1}/{n_perms}  "
                  f"null mean so far: {100 * np.mean(fracs):.2f}%")

    for f in [tmp_shuffled, tmp_closest]:
        if os.path.exists(f):
            os.remove(f)

    return np.array(fracs)


# ── Main ─────────────────────────────────────────────────────────────────────

def main(n_perms: int, seed: int) -> None:
    # Independent RNGs so each set of permutations is reproducible on its own
    rng_free = np.random.default_rng(seed)
    rng_sd   = np.random.default_rng(seed + 1)

    glyco_df    = pd.read_csv(GLYCO_TSV, sep="\t")
    glyco_names = set(glyco_df["hgnc_symbol"].dropna())

    obs_frac, obs_n_glyco, obs_n_total = glyco_fraction(CLOSEST_BED, glyco_names)
    print(f"Observed: {obs_n_glyco}/{obs_n_total} unique genes are glycoproteins "
          f"({100 * obs_frac:.2f}%)\n")

    with tempfile.TemporaryDirectory() as tmp:
        autosome_genome = make_autosome_genome(GENOME_FILE, tmp)
        sd_excl_bed     = prepare_sd_excl(SD_BED, tmp)

        merged_bed = os.path.join(tmp, "merged.bed")
        n_merged   = merge_bed(INPUT_BED, merged_bed)
        n_input    = sum(1 for _ in open(INPUT_BED))
        print(f"Input regions: {n_input} → {n_merged} after merge\n")

        print(f"Running {n_perms} unrestricted permutations…")
        null_free = run_permutations(
            merged_bed, autosome_genome, glyco_names, n_perms,
            rng_free, label="unrestricted",
        )

        print(f"\nRunning {n_perms} SD-excluded permutations…")
        null_sd = run_permutations(
            merged_bed, autosome_genome, glyco_names, n_perms,
            rng_sd, label="SD-excluded", excl_bed=sd_excl_bed,
        )

    # ── Report ────────────────────────────────────────────────────────────────

    for label, null in [("Unrestricted", null_free), ("SD-excluded", null_sd)]:
        p = np.mean(null >= obs_frac)
        p_str = f"{p:.4f}" if p > 0 else f"< {1/n_perms:.4f}"
        print(f"\n{label}:")
        print(f"  Null  mean={100*null.mean():.2f}%  sd={100*null.std():.2f}%")
        print(f"  Observed: {100*obs_frac:.2f}%")
        print(f"  Empirical p (>= observed): {p_str}  (N={n_perms})")

    # ── Save TSV ──────────────────────────────────────────────────────────────

    results_tsv = OUT_DIR / "glycoprotein_permutation_results.tsv"
    pd.DataFrame({
        "null_free_glyco_fraction": null_free,
        "null_sd_glyco_fraction":   null_sd,
    }).to_csv(results_tsv, sep="\t", index=False)
    print(f"\nResults saved to: {results_tsv}")

    # ── Plot ──────────────────────────────────────────────────────────────────

    plt.rcParams.update({
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5), sharey=True)

    configs = [
        (axes[0], null_free, "#7aabd9", "Unrestricted shuffle"),
        (axes[1], null_sd,   "#85c1a3", "SD-excluded shuffle"),
    ]

    for ax, null, color, title in configs:
        p = np.mean(null >= obs_frac)
        p_str = f"p = {p:.4f}" if p > 0 else f"p < {1/n_perms:.4f}"

        ax.hist(null * 100, bins=40, color=color, edgecolor="white",
                linewidth=0.5, label=f"Null ({n_perms} perms)")
        ax.axvline(obs_frac * 100, color="#e07070", lw=2,
                   label=f"Observed: {100*obs_frac:.1f}%")
        ax.axvline(null.mean() * 100, color="#444", lw=1.2, ls="--",
                   label=f"Null mean: {100*null.mean():.1f}%")

        ax.text(0.97, 0.95, p_str, transform=ax.transAxes,
                ha="right", va="top", fontsize=10,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="grey", alpha=0.8))

        ax.set_xlabel("Glycoprotein fraction of unique nearest genes (%)")
        ax.set_title(title)
        ax.legend(fontsize=9)

    axes[0].set_ylabel("Count (permutations)")
    fig.suptitle(
        "Permutation test: glycoprotein enrichment vs. genomic null\n"
        f"Observed: {obs_n_glyco}/{obs_n_total} genes = {100*obs_frac:.1f}%  |  "
        f"GENCODE background: 22.4%",
        fontsize=11,
    )
    fig.tight_layout()

    out_pdf = OUT_DIR / "glycoprotein_permutation_test.pdf"
    fig.savefig(out_pdf, bbox_inches="tight")
    print(f"Plot saved to: {out_pdf}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-perms", type=int, default=1000)
    parser.add_argument("--seed",    type=int, default=42)
    args = parser.parse_args()
    main(args.n_perms, args.seed)
