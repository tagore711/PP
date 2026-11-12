import os
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# CONFIGURATION
# ============================================================

INPUT_FILE = r"results\combined\all_chromosomes_region_counts.csv"

OUTPUT_DIR = r"results\visualizations\theta_below_10000"

LOW_SNP_LIMIT = 10000


# ============================================================
# CREATE OUTPUT DIRECTORY
# ============================================================

os.makedirs(OUTPUT_DIR, exist_ok=True)


# ============================================================
# READ DATA
# ============================================================

print("=" * 70)
print("THETA REGION SNP ANALYSIS")
print("=" * 70)

print(f"\nReading: {INPUT_FILE}")

df = pd.read_csv(INPUT_FILE)

print(f"Total rows: {len(df)}")
print(f"Columns: {list(df.columns)}")


# ============================================================
# CHECK REQUIRED COLUMNS
# ============================================================

required_columns = [
    "chromosome",
    "region_number",
    "snp_count"
]

missing = [
    col for col in required_columns
    if col not in df.columns
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
# CHROMOSOME SORTING
# ============================================================

def chromosome_number(chromosome):
    return int(
        str(chromosome)
        .replace("chr", "")
    )


chromosomes = sorted(
    df["chromosome"].unique(),
    key=chromosome_number
)


# ============================================================
# FIND REGIONS BELOW 10,000 SNPs
# ============================================================

low_snp_df = df[
    df["snp_count"] < LOW_SNP_LIMIT
].copy()

low_snp_df = low_snp_df.sort_values(
    ["chromosome", "region_number"],
    key=lambda col:
        col.map(chromosome_number)
        if col.name == "chromosome"
        else col
)


# ============================================================
# SAVE LOW-SNP REGIONS
# ============================================================

low_snp_file = os.path.join(
    OUTPUT_DIR,
    "theta_regions_snp_less_than_10000.csv"
)

low_snp_df.to_csv(
    low_snp_file,
    index=False
)


# ============================================================
# PRINT SUMMARY
# ============================================================

print("\n" + "=" * 70)
print("REGIONS WITH SNP COUNT < 10,000")
print("=" * 70)

print(
    f"\nTotal regions below {LOW_SNP_LIMIT}: "
    f"{len(low_snp_df)}"
)

print(f"Saved: {low_snp_file}")


# ============================================================
# CREATE 22 CHROMOSOME GRAPHS
# ============================================================

summary = []


for chromosome in chromosomes:

    # Get this chromosome
    chrom_df = df[
        df["chromosome"] == chromosome
    ].copy()

    # Sort θ regions
    chrom_df = chrom_df.sort_values(
        "region_number"
    )

    # Regions below 10,000
    low_df = chrom_df[
        chrom_df["snp_count"] < LOW_SNP_LIMIT
    ]

    # --------------------------------------------------------
    # CREATE FIGURE
    # --------------------------------------------------------

    plt.figure(
        figsize=(14, 7)
    )

    # Every θ region
    plt.plot(
        chrom_df["region_number"],
        chrom_df["snp_count"],
        marker="o",
        markersize=3,
        linewidth=1,
        label="SNP count"
    )

    # Highlight < 10,000
    if len(low_df) > 0:

        plt.scatter(
            low_df["region_number"],
            low_df["snp_count"],
            s=55,
            zorder=5,
            label="SNP count < 10,000"
        )

    # 10,000 threshold
    plt.axhline(
        LOW_SNP_LIMIT,
        linestyle="--",
        linewidth=1.5,
        label="10,000 threshold"
    )

    # --------------------------------------------------------
    # LABELS
    # --------------------------------------------------------

    plt.title(
        f"{chromosome}: SNP Count Across θ Regions",
        fontsize=16,
        fontweight="bold"
    )

    plt.xlabel(
        "θ Region Number",
        fontsize=12
    )

    plt.ylabel(
        "SNP Count",
        fontsize=12
    )

    plt.grid(
        True,
        alpha=0.3
    )

    plt.legend()

    # --------------------------------------------------------
    # X AXIS
    # --------------------------------------------------------

    number_of_regions = len(chrom_df)

    if number_of_regions <= 100:

        plt.xticks(
            chrom_df["region_number"],
            rotation=90,
            fontsize=7
        )

    # --------------------------------------------------------
    # SAVE
    # --------------------------------------------------------

    output_file = os.path.join(
        OUTPUT_DIR,
        f"{chromosome}_theta_snp_count.png"
    )

    plt.tight_layout()

    plt.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

    # --------------------------------------------------------
    # SUMMARY
    # --------------------------------------------------------

    summary.append({
        "chromosome": chromosome,
        "total_theta_regions": len(chrom_df),
        "regions_below_10000": len(low_df),
        "minimum_snp_count": int(
            chrom_df["snp_count"].min()
        ),
        "maximum_snp_count": int(
            chrom_df["snp_count"].max()
        )
    })

    print(
        f"{chromosome}: "
        f"{len(chrom_df)} θ regions | "
        f"{len(low_df)} below 10,000"
    )


# ============================================================
# SAVE SUMMARY
# ============================================================

summary_df = pd.DataFrame(summary)

summary_file = os.path.join(
    OUTPUT_DIR,
    "theta_less_than_10000_summary.csv"
)

summary_df.to_csv(
    summary_file,
    index=False
)


# ============================================================
# FINAL MESSAGE
# ============================================================

print("\n" + "=" * 70)
print("COMPLETED SUCCESSFULLY")
print("=" * 70)

print(
    f"\n22 chromosome graphs saved to:"
    f"\n{OUTPUT_DIR}"
)

print(
    f"\nLow-SNP θ regions:"
    f"\n{low_snp_file}"
)

print(
    f"\nSummary:"
    f"\n{summary_file}"
)

print("\nDone.")