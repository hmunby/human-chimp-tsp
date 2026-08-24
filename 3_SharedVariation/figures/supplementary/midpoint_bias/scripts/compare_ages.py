import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

GENERATION_TIME = 29  # years per generation
INFILE = snakemake.input[0]
OUTFILE = snakemake.output[0]
N_BINS = 50  # log-spaced age bins

with open(INFILE, "rb") as f:
    data = pickle.load(f)

actual_ages, midpoint_ages = zip(*data)
actual_years = np.array(actual_ages) * GENERATION_TIME
midpoint_years = np.array(midpoint_ages) * GENERATION_TIME

# Log-spaced bins across the range of actual ages
bin_edges = np.logspace(
    np.log10(actual_years.min()),
    np.log10(actual_years.max()),
    N_BINS + 1
)
bin_indices = np.digitize(actual_years, bin_edges)

bin_actual_mean = []
bin_midpoint_mean = []
bin_midpoint_median = []
bin_midpoint_q25 = []
bin_midpoint_q75 = []
bin_centers = []

for i in range(1, N_BINS + 1):
    mask = bin_indices == i
    if mask.sum() < 10:
        continue
    mp = midpoint_years[mask]
    bin_actual_mean.append(actual_years[mask].mean())
    bin_midpoint_mean.append(mp.mean())
    bin_midpoint_median.append(np.median(mp))
    bin_midpoint_q25.append(np.percentile(mp, 25))
    bin_midpoint_q75.append(np.percentile(mp, 75))
    bin_centers.append(actual_years[mask].mean())

bin_actual_mean = np.array(bin_actual_mean)
bin_midpoint_mean = np.array(bin_midpoint_mean)
bin_midpoint_median = np.array(bin_midpoint_median)
bin_midpoint_q25 = np.array(bin_midpoint_q25)
bin_midpoint_q75 = np.array(bin_midpoint_q75)
bin_centers = np.array(bin_centers)

fig, axes = plt.subplots(2, 1, figsize=(8, 9))

# Panel 1: actual vs midpoint age
ax = axes[0]
ax.fill_between(bin_centers, bin_midpoint_q25, bin_midpoint_q75,
                color="steelblue", alpha=0.2, label="IQR")
ax.plot(bin_actual_mean, bin_actual_mean, "k--", lw=1, label="1:1 line")
ax.plot(bin_centers, bin_midpoint_mean, color="steelblue", lw=2, label="Mean midpoint age")
ax.plot(bin_centers, bin_midpoint_median, color="steelblue", lw=2, ls=":", label="Median midpoint age")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Actual mutation age (years)")
ax.set_ylabel("Midpoint branch age (years)")
ax.set_title("Actual vs. midpoint branch age")
ax.legend()
ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:,.0f}"))
ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:,.0f}"))

# Panel 2: ratio of midpoint to actual age
ax = axes[1]
ratio_mean = bin_midpoint_mean / bin_actual_mean
ratio_median = bin_midpoint_median / bin_actual_mean
ratio_q25 = bin_midpoint_q25 / bin_actual_mean
ratio_q75 = bin_midpoint_q75 / bin_actual_mean
ax.fill_between(bin_centers, ratio_q25, ratio_q75,
                color="steelblue", alpha=0.2, label="IQR")
ax.axhline(1, color="k", lw=1, ls="--")
ax.plot(bin_centers, ratio_mean, color="steelblue", lw=2, label="Mean midpoint / actual")
ax.plot(bin_centers, ratio_median, color="steelblue", lw=2, ls=":", label="Median midpoint / actual")
ax.set_xscale("log")
ax.set_xlabel("Actual mutation age (years)")
ax.set_ylabel("Midpoint age / actual age")
ax.set_title("Midpoint age as a fraction of actual age")
ax.legend()
ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:,.0f}"))

plt.tight_layout()
plt.savefig(OUTFILE, dpi=150)
