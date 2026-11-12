
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# --------------------------------------------------
# PATHS
# --------------------------------------------------
ROOT = Path(__file__).resolve().parent.parent

INPUT_FILE = (
    ROOT / "results" / "combined"
    / "all_chromosomes_region_counts.csv"
)

OUTPUT_DIR = ROOT / "results" / "visualizations"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------
# LOAD DATA
# --------------------------------------------------
df = pd.read_csv(INPUT_FILE)

if "snp_count" not in df.columns:
    raise ValueError(
        f"snp_count column missing. Available: {list(df.columns)}"
    )

counts = pd.to_numeric(
    df["snp_count"], errors="raise"
).dropna()

if counts.empty:
    raise ValueError("No SNP counts found.")

mean_count = counts.mean()
median_count = counts.median()

percentiles = {
    "P1": counts.quantile(0.01),
    "P5": counts.quantile(0.05),
    "P10": counts.quantile(0.10),
    "P25": counts.quantile(0.25),
    "P50 (Median)": counts.quantile(0.50),
}

# --------------------------------------------------
# GRAPH 1: HISTOGRAM
# --------------------------------------------------
fig, ax = plt.subplots(figsize=(12, 6))

ax.hist(
    counts,
    bins=50,
    edgecolor="black",
    alpha=0.8
)

ax.axvline(
    mean_count,
    linestyle="--",
    linewidth=2,
    label=f"Mean = {mean_count:,.0f}"
)

ax.axvline(
    median_count,
    linestyle=":",
    linewidth=2,
    label=f"Median = {median_count:,.0f}"
)

for label in ["P5", "P10"]:
    value = percentiles[label]
    ax.axvline(
        value,
        linestyle="--",
        linewidth=1.5,
        label=f"{label} = {value:,.0f}"
    )

ax.set_title("Distribution of SNP Counts Across θ Regions")
ax.set_xlabel("SNP Count per θ Region (1 Mbp)")
ax.set_ylabel("Number of θ Regions")
ax.legend()

fig.tight_layout()
fig.savefig(
    OUTPUT_DIR / "theta_snp_histogram.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close(fig)

# --------------------------------------------------
# GRAPH 2: SORTED SNP COUNTS
# --------------------------------------------------
sorted_counts = np.sort(counts.to_numpy())
region_ranks = np.arange(1, len(sorted_counts) + 1)

fig, ax = plt.subplots(figsize=(13, 6))

ax.plot(
    region_ranks,
    sorted_counts,
    linewidth=1.5
)

ax.axhline(
    mean_count,
    linestyle="--",
    linewidth=2,
    label=f"Mean = {mean_count:,.0f}"
)

for label in ["P5", "P10", "P25", "P50 (Median)"]:
    value = percentiles[label]
    ax.axhline(
        value,
        linestyle=":",
        linewidth=1.5,
        label=f"{label} = {value:,.0f}"
    )

ax.set_title("Sorted SNP Counts Across All θ Regions")
ax.set_xlabel("θ Regions (Sorted from Lowest to Highest SNP Count)")
ax.set_ylabel("SNP Count per θ Region")
ax.legend()

fig.tight_layout()
fig.savefig(
    OUTPUT_DIR / "theta_snp_sorted_curve.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close(fig)

# --------------------------------------------------
# PRINT CANDIDATE THRESHOLDS
# --------------------------------------------------
print("\n===== θ REGION SNP THRESHOLD ANALYSIS =====")
print(f"Total regions: {len(counts):,}")
print(f"Mean SNP count: {mean_count:,.2f}")
print(f"Median SNP count: {median_count:,.2f}")
print(f"Minimum SNP count: {counts.min():,.0f}")
print(f"Maximum SNP count: {counts.max():,.0f}")

print("\nCandidate thresholds:")
for label, value in percentiles.items():
    retained = (counts >= value).sum()
    excluded = (counts < value).sum()

    print(
        f"{label}: {value:,.0f} SNPs | "
        f"Retained: {retained:,} | "
        f"Below threshold: {excluded:,}"
    )

print(f"\nGraphs saved to: {OUTPUT_DIR}")
