from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch


# ============================================================
# CONFIGURATION
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

INPUT_FILE = (
    PROJECT_ROOT
    / "results"
    / "combined"
    / "all_chromosomes_region_counts.csv"
)

OUTPUT_DIR = (
    PROJECT_ROOT
    / "results"
    / "visualizations"
    / "theta_analysis"
)

STATS_DIR = (
    PROJECT_ROOT
    / "results"
    / "statistics"
)

THETA_SIZE = 1_000_000

# IMPORTANT:
# "fewer than 1000 SNPs" means STRICTLY < 1000.
# Therefore:
#   0-999   -> REMOVE CANDIDATE
#   1000+   -> KEEP
LOW_THRESHOLD = 1000


# ============================================================
# SETUP
# ============================================================

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
STATS_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# LOAD DATA
# ============================================================

df = pd.read_csv(INPUT_FILE)

required_columns = [
    "chromosome",
    "chromosome_length_bp",
    "length_source",
    "theta_size_bp",
    "number_of_theta_regions",
    "region_number",
    "region_start_bp",
    "region_end_bp",
    "snp_count",
]

missing = [c for c in required_columns if c not in df.columns]

if missing:
    raise ValueError(
        f"Missing required columns: {missing}"
    )


# ============================================================
# VALIDATION
# ============================================================

if not (df["theta_size_bp"] == THETA_SIZE).all():
    raise ValueError(
        "Input contains theta sizes other than 1,000,000 bp."
    )

chromosomes = [f"chr{i}" for i in range(1, 23)]

df["chromosome"] = pd.Categorical(
    df["chromosome"],
    categories=chromosomes,
    ordered=True,
)

df = df.sort_values(
    ["chromosome", "region_number"]
).reset_index(drop=True)


# ============================================================
# CLASSIFICATION
# ============================================================

# STRICT interpretation:
#
# < 1000 SNPs = removal candidate
# >= 1000 SNPs = keep

df["filter_category"] = np.where(
    df["snp_count"] < LOW_THRESHOLD,
    "REMOVE_CANDIDATE",
    "KEEP",
)


# ============================================================
# CREATE EXACT REMOVAL-CANDIDATE CSV
# ============================================================

remove_df = df[
    df["snp_count"] < LOW_THRESHOLD
].copy()

remove_df = remove_df[
    [
        "chromosome",
        "region_number",
        "region_start_bp",
        "region_end_bp",
        "snp_count",
        "theta_size_bp",
        "filter_category",
    ]
]

remove_df.to_csv(
    STATS_DIR / "regions_below_1000_snps.csv",
    index=False,
)


# ============================================================
# CREATE CHROMOSOME-WISE SUMMARY
# ============================================================

summary_rows = []

for chrom in chromosomes:

    chrom_df = df[
        df["chromosome"] == chrom
    ]

    total_regions = len(chrom_df)

    remove_count = (
        chrom_df["snp_count"] < LOW_THRESHOLD
    ).sum()

    zero_count = (
        chrom_df["snp_count"] == 0
    ).sum()

    one_to_999_count = (
        (chrom_df["snp_count"] >= 1)
        & (chrom_df["snp_count"] < LOW_THRESHOLD)
    ).sum()

    keep_count = (
        chrom_df["snp_count"] >= LOW_THRESHOLD
    ).sum()

    summary_rows.append(
        {
            "chromosome": chrom,
            "total_existing_regions": total_regions,
            "regions_below_1000": remove_count,
            "regions_with_0_snps": zero_count,
            "regions_with_1_to_999_snps": one_to_999_count,
            "regions_1000_or_more": keep_count,
            "percent_below_1000": (
                remove_count / total_regions * 100
                if total_regions > 0
                else 0
            ),
        }
    )


summary_df = pd.DataFrame(summary_rows)

summary_df.to_csv(
    STATS_DIR / "regions_below_1000_summary_by_chromosome.csv",
    index=False,
)


# ============================================================
# OVERALL SUMMARY
# ============================================================

total_existing = len(df)

total_remove = (
    df["snp_count"] < LOW_THRESHOLD
).sum()

total_zero = (
    df["snp_count"] == 0
).sum()

total_1_to_999 = (
    (df["snp_count"] >= 1)
    & (df["snp_count"] < LOW_THRESHOLD)
).sum()

total_keep = (
    df["snp_count"] >= LOW_THRESHOLD
).sum()

overall_summary = pd.DataFrame(
    [
        {
            "total_existing_theta_regions": total_existing,
            "regions_below_1000": total_remove,
            "regions_with_0_snps": total_zero,
            "regions_with_1_to_999_snps": total_1_to_999,
            "regions_1000_or_more": total_keep,
            "percent_below_1000": (
                total_remove / total_existing * 100
            ),
            "percent_1000_or_more": (
                total_keep / total_existing * 100
            ),
        }
    ]
)

overall_summary.to_csv(
    STATS_DIR / "regions_below_1000_overall_summary.csv",
    index=False,
)


# ============================================================
# PRINT SUMMARY
# ============================================================

print()
print("=" * 70)
print("THETA REGION FILTERING ANALYSIS")
print("=" * 70)

print(f"Theta size              : {THETA_SIZE:,} bp")
print(f"Removal threshold       : < {LOW_THRESHOLD:,} SNPs")
print()
print(f"Existing theta regions  : {total_existing:,}")
print(f"< 1000 SNPs             : {total_remove:,}")
print(f"    Exactly 0           : {total_zero:,}")
print(f"    1-999               : {total_1_to_999:,}")
print(f">= 1000 SNPs            : {total_keep:,}")
print()
print(
    f"Removal candidates      : "
    f"{total_remove / total_existing * 100:.2f}%"
)
print(
    f"Regions to keep         : "
    f"{total_keep / total_existing * 100:.2f}%"
)

print()
print("Output files:")
print(
    "  regions_below_1000_snps.csv"
)
print(
    "  regions_below_1000_summary_by_chromosome.csv"
)
print(
    "  regions_below_1000_overall_summary.csv"
)


# ============================================================
# HEATMAP FUNCTION
# ============================================================

def create_filter_heatmap(
    start_region,
    end_region,
    output_filename,
):

    region_numbers = list(
        range(start_region, end_region + 1)
    )

    matrix = np.full(
        (len(region_numbers), len(chromosomes)),
        np.nan,
    )

    # 0 = below 1000 -> REMOVE
    # 1 = >= 1000   -> KEEP

    for col, chrom in enumerate(chromosomes):

        chrom_df = df[
            df["chromosome"] == chrom
        ]

        lookup = dict(
            zip(
                chrom_df["region_number"],
                chrom_df["snp_count"],
            )
        )

        for row, region in enumerate(region_numbers):

            if region in lookup:

                snp_count = lookup[region]

                if snp_count < LOW_THRESHOLD:
                    matrix[row, col] = 0
                else:
                    matrix[row, col] = 1

            # Non-existent regions remain NaN.
            # They will be displayed as white.

    # ========================================================
    # COLORS
    # ========================================================

    cmap = ListedColormap(
        [
            "#D73027",  # RED = remove candidate
            "#4575B4",  # BLUE = keep
        ]
    )

    cmap.set_bad("#FFFFFF")

    masked_matrix = np.ma.masked_invalid(matrix)

    # ========================================================
    # FIGURE
    # ========================================================

    fig_height = max(
        9,
        len(region_numbers) * 0.085
    )

    fig, ax = plt.subplots(
        figsize=(15, fig_height),
        dpi=300,
    )

    ax.imshow(
        masked_matrix,
        aspect="auto",
        interpolation="nearest",
        cmap=cmap,
        vmin=0,
        vmax=1,
    )

    # ========================================================
    # AXES
    # ========================================================

    ax.set_xticks(
        np.arange(len(chromosomes))
    )

    ax.set_xticklabels(
        chromosomes,
        fontsize=10,
        fontweight="bold",
    )

    ax.set_xlabel(
        "Chromosome",
        fontsize=12,
        fontweight="bold",
    )

    # Only show every 10th region label to avoid congestion.
    tick_positions = []

    for i, region in enumerate(region_numbers):

        if (
            region % 10 == 0
            or region == start_region
            or region == end_region
        ):
            tick_positions.append(i)

    ax.set_yticks(tick_positions)

    ax.set_yticklabels(
        [
            str(region_numbers[i])
            for i in tick_positions
        ],
        fontsize=9,
    )

    ax.set_ylabel(
        "Theta region number",
        fontsize=12,
        fontweight="bold",
    )

    # ========================================================
    # GRID
    # ========================================================

    ax.set_xticks(
        np.arange(-0.5, len(chromosomes), 1),
        minor=True,
    )

    ax.grid(
        which="minor",
        axis="x",
        linewidth=0.5,
        alpha=0.35,
    )

    ax.tick_params(
        which="minor",
        bottom=False,
        left=False,
    )

    # ========================================================
    # TITLE
    # ========================================================

    ax.set_title(
        f"Theta Regions {start_region}–{end_region}: "
        f"Regions Below 1,000 SNP Records",
        fontsize=16,
        fontweight="bold",
        pad=15,
    )

    # ========================================================
    # LEGEND
    # ========================================================

    legend_elements = [
        Patch(
            facecolor="#D73027",
            edgecolor="black",
            label="< 1,000 SNPs — REMOVE CANDIDATE",
        ),
        Patch(
            facecolor="#4575B4",
            edgecolor="black",
            label="≥ 1,000 SNPs — KEEP",
        ),
        Patch(
            facecolor="#FFFFFF",
            edgecolor="black",
            label="Region does not exist",
        ),
    ]

    ax.legend(
        handles=legend_elements,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.055),
        ncol=3,
        frameon=True,
        fontsize=10,
    )

    plt.tight_layout()

    output_path = (
        OUTPUT_DIR / output_filename
    )

    fig.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
    )

    plt.close(fig)

    print(
        f"Created: {output_path}"
    )


# ============================================================
# CREATE FINAL TWO HEATMAPS
# ============================================================

create_filter_heatmap(
    1,
    124,
    "theta_regions_1_124_filter_heatmap.png",
)

create_filter_heatmap(
    125,
    249,
    "theta_regions_125_249_filter_heatmap.png",
)


# ============================================================
# FINAL VALIDATION
# ============================================================

print()
print("=" * 70)
print("VALIDATION")
print("=" * 70)

print(
    f"Source records            : {len(df):,}"
)

print(
    f"Removal candidates        : {len(remove_df):,}"
)

print(
    f"Keep regions              : {total_keep:,}"
)

print(
    f"0-SNP regions             : {total_zero:,}"
)

print(
    f"1-999 SNP regions         : {total_1_to_999:,}"
)

print(
    f"1000+ SNP regions         : {total_keep:,}"
)

assert (
    total_remove
    + total_keep
    == total_existing
)

assert (
    total_zero
    + total_1_to_999
    == total_remove
)

assert (
    (df["snp_count"] < 1000)
    .sum()
    == len(remove_df)
)

print()
print("ALL VALIDATIONS PASSED.")
print("=" * 70)