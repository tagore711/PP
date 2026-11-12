import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


# ============================================================
# CONFIGURATION
# ============================================================

INPUT_FILE = r"results\combined\all_chromosomes_region_counts.csv"

OUTPUT_DIR = r"results\visualizations\theta_below_10000"

THRESHOLD = 10000


# ============================================================
# CREATE OUTPUT DIRECTORY
# ============================================================

os.makedirs(OUTPUT_DIR, exist_ok=True)


# ============================================================
# READ DATA
# ============================================================

print("=" * 75)
print("θ REGION IDENTIFICATION: SNP COUNT < 10,000")
print("=" * 75)

print(f"\nReading:")
print(INPUT_FILE)

df = pd.read_csv(INPUT_FILE)

print(f"\nTotal θ-region records: {len(df)}")

print("\nColumns:")
print(list(df.columns))


# ============================================================
# CHECK REQUIRED COLUMNS
# ============================================================

required_columns = [
    "chromosome",
    "region_number",
    "snp_count"
]

missing = [
    column
    for column in required_columns
    if column not in df.columns
]

if missing:
    raise ValueError(
        f"Missing required columns: {missing}"
    )


# ============================================================
# CLEAN DATA
# ============================================================

df["chromosome"] = df["chromosome"].astype(str)

df["region_number"] = pd.to_numeric(
    df["region_number"],
    errors="coerce"
)

df["snp_count"] = pd.to_numeric(
    df["snp_count"],
    errors="coerce"
)

df = df.dropna(
    subset=[
        "chromosome",
        "region_number",
        "snp_count"
    ]
)


# ============================================================
# CHROMOSOME ORDER
# ============================================================

chromosomes = [
    f"chr{i}"
    for i in range(1, 23)
]

# Keep only chromosomes actually present
chromosomes = [
    chromosome
    for chromosome in chromosomes
    if chromosome in df["chromosome"].unique()
]


# ============================================================
# IDENTIFY LOW-SNP θ REGIONS
# ============================================================

low_df = df[
    df["snp_count"] < THRESHOLD
].copy()

low_df = low_df.sort_values(
    ["chromosome", "region_number"]
)


print("\n" + "=" * 75)
print(f"TOTAL θ REGIONS WITH SNP COUNT < {THRESHOLD}")
print("=" * 75)

print(
    f"\nTotal: {len(low_df)}"
)


# ============================================================
# PRINT EXACT θ REGIONS
# ============================================================

print("\n" + "=" * 75)
print("EXACT θ REGIONS")
print("=" * 75)

for chromosome in chromosomes:

    chrom_low = low_df[
        low_df["chromosome"] == chromosome
    ]

    regions = (
        chrom_low["region_number"]
        .astype(int)
        .tolist()
    )

    print(
        f"\n{chromosome}: "
        f"{len(regions)} regions"
    )

    if regions:
        print(
            "θ = "
            + ", ".join(
                str(region)
                for region in regions
            )
        )
    else:
        print("θ = None")


# ============================================================
# SAVE EXACT LOW-SNP REGIONS
# ============================================================

exact_regions_file = os.path.join(
    OUTPUT_DIR,
    "theta_regions_below_10000_exact.csv"
)

low_df.to_csv(
    exact_regions_file,
    index=False
)

print("\nSaved exact region data:")
print(exact_regions_file)


# ============================================================
# CREATE SIMPLE CHROMOSOME → θ REGION TABLE
# ============================================================

rows = []

for chromosome in chromosomes:

    chrom_low = low_df[
        low_df["chromosome"] == chromosome
    ]

    regions = (
        chrom_low["region_number"]
        .astype(int)
        .tolist()
    )

    rows.append({
        "chromosome": chromosome,
        "number_of_regions_below_10000": len(regions),
        "theta_regions_below_10000": ", ".join(
            str(region)
            for region in regions
        ) if regions else "None"
    })


region_table = pd.DataFrame(rows)

region_table_file = os.path.join(
    OUTPUT_DIR,
    "theta_regions_below_10000_by_chromosome.csv"
)

region_table.to_csv(
    region_table_file,
    index=False
)

print("\nSaved chromosome-wise θ-region table:")
print(region_table_file)


# ============================================================
# CREATE θ REGION MATRIX
# ============================================================

max_theta = int(
    df["region_number"].max()
)

matrix = np.full(
    (max_theta, len(chromosomes)),
    np.nan
)


for x, chromosome in enumerate(chromosomes):

    chrom_df = df[
        df["chromosome"] == chromosome
    ]

    for _, row in chrom_df.iterrows():

        theta = int(
            row["region_number"]
        )

        snp_count = row["snp_count"]

        if snp_count < THRESHOLD:

            matrix[
                theta - 1,
                x
            ] = 1

        else:

            matrix[
                theta - 1,
                x
            ] = 0


# ============================================================
# PLOT θ REGION HEATMAP
# ============================================================

print("\nCreating θ-region heatmap...")


fig, ax = plt.subplots(
    figsize=(18, 15)
)


# Mask NaN values so nonexistent θ regions are blank
masked_matrix = np.ma.masked_invalid(
    matrix
)


image = ax.imshow(
    masked_matrix,
    aspect="auto",
    interpolation="none",
    origin="lower"
)


# ------------------------------------------------------------
# X AXIS
# ------------------------------------------------------------

ax.set_xticks(
    np.arange(len(chromosomes))
)

ax.set_xticklabels(
    chromosomes,
    fontsize=11
)


ax.set_xlabel(
    "Chromosome",
    fontsize=13,
    fontweight="bold"
)


# ------------------------------------------------------------
# Y AXIS
# ------------------------------------------------------------

ax.set_ylabel(
    "θ Region Number",
    fontsize=13,
    fontweight="bold"
)


# Show θ numbers
# For very large ranges, use sensible spacing

if max_theta <= 100:

    yticks = np.arange(
        1,
        max_theta + 1
    )

else:

    step = 10

    yticks = np.arange(
        1,
        max_theta + 1,
        step
    )


ax.set_yticks(
    yticks - 1
)

ax.set_yticklabels(
    yticks,
    fontsize=9
)


# ------------------------------------------------------------
# TITLE
# ------------------------------------------------------------

ax.set_title(
    "θ Regions with SNP Count < 10,000 Across Chromosomes",
    fontsize=17,
    fontweight="bold",
    pad=15
)


# ------------------------------------------------------------
# GRID
# ------------------------------------------------------------

ax.set_xticks(
    np.arange(
        -0.5,
        len(chromosomes),
        1
    ),
    minor=True
)

ax.set_yticks(
    np.arange(
        -0.5,
        max_theta,
        1
    ),
    minor=True
)

ax.grid(
    which="minor",
    linewidth=0.15,
    alpha=0.3
)

ax.tick_params(
    which="minor",
    bottom=False,
    left=False
)


# ------------------------------------------------------------
# COLORBAR
# ------------------------------------------------------------

colorbar = plt.colorbar(
    image,
    ax=ax,
    ticks=[0, 1],
    fraction=0.025,
    pad=0.02
)

colorbar.ax.set_yticklabels(
    [
        f">= {THRESHOLD:,} SNPs",
        f"< {THRESHOLD:,} SNPs"
    ]
)

colorbar.set_label(
    "SNP Count Category",
    fontsize=11
)


# ------------------------------------------------------------
# LAYOUT
# ------------------------------------------------------------

plt.tight_layout()


# ============================================================
# SAVE HEATMAP
# ============================================================

heatmap_file = os.path.join(
    OUTPUT_DIR,
    "theta_regions_below_10000_heatmap.png"
)

plt.savefig(
    heatmap_file,
    dpi=300,
    bbox_inches="tight"
)

plt.close()


# ============================================================
# FINAL SUMMARY
# ============================================================

print("\n" + "=" * 75)
print("COMPLETED")
print("=" * 75)

print("\nHeatmap:")
print(heatmap_file)

print("\nExact θ-region CSV:")
print(exact_regions_file)

print("\nChromosome-wise θ-region table:")
print(region_table_file)

print(
    f"\nTotal θ regions < {THRESHOLD:,}: "
    f"{len(low_df)}"
)

print("\nDone.")