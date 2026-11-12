from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


# ============================================================
# PATHS
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

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
STATS_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# CONFIGURATION
# ============================================================

CHROMOSOME = "chr1"
THETA_SIZE = 1_000_000
THRESHOLD = 1_000


# ============================================================
# LOAD DATA
# ============================================================

df = pd.read_csv(INPUT_FILE)

chr1 = df[
    df["chromosome"] == CHROMOSOME
].copy()

chr1 = chr1.sort_values(
    "region_number"
).reset_index(drop=True)


# ============================================================
# VALIDATION
# ============================================================

expected_regions = list(range(1, 250))

actual_regions = (
    chr1["region_number"]
    .astype(int)
    .tolist()
)

if actual_regions != expected_regions:
    raise ValueError(
        "Chromosome 1 does not contain theta regions 1–249."
    )

if not (
    chr1["theta_size_bp"] == THETA_SIZE
).all():
    raise ValueError(
        "Theta size is not 1 Mbp."
    )


# ============================================================
# CLASSIFICATION
# ============================================================

chr1["category"] = chr1["snp_count"].apply(
    lambda x:
        "Remove Candidate"
        if x < THRESHOLD
        else "Keep"
)


remove_candidates = chr1[
    chr1["snp_count"] < THRESHOLD
].copy()

keep_regions = chr1[
    chr1["snp_count"] >= THRESHOLD
].copy()


# ============================================================
# SAVE EXACT CANDIDATE LIST
# ============================================================

candidate_output = remove_candidates[
    [
        "chromosome",
        "region_number",
        "region_start_bp",
        "region_end_bp",
        "snp_count",
        "category",
    ]
].copy()

candidate_output.to_csv(
    STATS_DIR / "chr1_regions_below_1000_snps.csv",
    index=False,
)


# ============================================================
# PRINT SUMMARY
# ============================================================

print()
print("=" * 75)
print("CHROMOSOME 1 — THETA REGION FILTERING")
print("=" * 75)

print(f"Chromosome             : {CHROMOSOME}")
print(f"Theta size             : {THETA_SIZE:,} bp")
print(f"Total theta regions    : {len(chr1)}")
print(
    f"Regions < 1,000 SNPs  : "
    f"{len(remove_candidates)}"
)
print(
    f"Regions >= 1,000 SNPs : "
    f"{len(keep_regions)}"
)

print()
print("REMOVE CANDIDATES")
print("-" * 75)

if len(remove_candidates) == 0:
    print("No theta regions contain fewer than 1,000 SNPs.")
else:
    print(
        remove_candidates[
            [
                "region_number",
                "region_start_bp",
                "region_end_bp",
                "snp_count",
            ]
        ].to_string(index=False)
    )


# ============================================================
# COLORS
# ============================================================

REMOVE_COLOR = "#D73027"
KEEP_COLOR = "#4575B4"
THRESHOLD_COLOR = "#222222"


# ============================================================
# GRAPH 1 — CLEAN OVERVIEW
# ============================================================

fig, ax = plt.subplots(
    figsize=(20, 8),
    dpi=300,
)

x = chr1["region_number"]
y = chr1["snp_count"]

colors = [
    REMOVE_COLOR if value < THRESHOLD else KEEP_COLOR
    for value in y
]

ax.bar(
    x,
    y,
    color=colors,
    width=0.8,
    edgecolor="none",
)


# Threshold line
ax.axhline(
    THRESHOLD,
    color=THRESHOLD_COLOR,
    linestyle="--",
    linewidth=2,
    label="1,000 SNP threshold",
)


# Threshold annotation
ax.text(
    249,
    THRESHOLD + 150,
    "1,000 SNP threshold",
    ha="right",
    va="bottom",
    fontsize=10,
    fontweight="bold",
)


# ------------------------------------------------------------
# X AXIS — reduce clutter
# ------------------------------------------------------------

# Major labels every 10 regions
major_ticks = list(range(10, 250, 10))

if 249 not in major_ticks:
    major_ticks.append(249)

ax.set_xticks(major_ticks)

ax.set_xticklabels(
    [str(x) for x in major_ticks],
    fontsize=9,
)

ax.set_xlim(0, 250)

ax.set_xlabel(
    "Theta region number (1 Mbp per region)",
    fontsize=12,
    fontweight="bold",
)

ax.set_ylabel(
    "SNP record count",
    fontsize=12,
    fontweight="bold",
)


# ------------------------------------------------------------
# TITLE
# ------------------------------------------------------------

ax.set_title(
    "Chromosome 1 — SNP Records Across Theta Regions",
    fontsize=17,
    fontweight="bold",
    pad=15,
)

ax.grid(
    axis="y",
    linestyle="--",
    alpha=0.25,
)

ax.set_axisbelow(True)


# ------------------------------------------------------------
# LEGEND
# ------------------------------------------------------------

legend_elements = [
    Patch(
        facecolor=REMOVE_COLOR,
        label="< 1,000 SNPs — Remove Candidate",
    ),
    Patch(
        facecolor=KEEP_COLOR,
        label="≥ 1,000 SNPs — Keep",
    ),
]

ax.legend(
    handles=legend_elements,
    loc="upper left",
    fontsize=10,
    frameon=True,
)


plt.tight_layout()

overview_path = (
    OUTPUT_DIR
    / "chr1_theta_snp_overview.png"
)

plt.savefig(
    overview_path,
    dpi=300,
    bbox_inches="tight",
    facecolor="white",
)

plt.close()


# ============================================================
# GRAPH 2 — CANDIDATE-ONLY ZOOM
# ============================================================

if len(remove_candidates) > 0:

    fig, ax = plt.subplots(
        figsize=(14, 7),
        dpi=300,
    )

    candidate_x = (
        remove_candidates["region_number"]
        .astype(int)
    )

    candidate_y = (
        remove_candidates["snp_count"]
        .astype(int)
    )

    bars = ax.bar(
        range(len(remove_candidates)),
        candidate_y,
        color=REMOVE_COLOR,
        width=0.65,
        edgecolor="black",
        linewidth=0.5,
    )


    # --------------------------------------------------------
    # X LABELS
    # --------------------------------------------------------

    labels = [
        f"θ{region}"
        for region in candidate_x
    ]

    ax.set_xticks(
        range(len(labels))
    )

    ax.set_xticklabels(
        labels,
        fontsize=10,
        fontweight="bold",
    )

    ax.set_xlabel(
        "Theta region number",
        fontsize=12,
        fontweight="bold",
    )


    # --------------------------------------------------------
    # Y AXIS
    # --------------------------------------------------------

    ax.set_ylabel(
        "SNP record count",
        fontsize=12,
        fontweight="bold",
    )

    ax.set_ylim(
        0,
        max(
            THRESHOLD * 1.05,
            max(candidate_y) * 1.25,
        ),
    )


    # --------------------------------------------------------
    # THRESHOLD
    # --------------------------------------------------------

    ax.axhline(
        THRESHOLD,
        color=THRESHOLD_COLOR,
        linestyle="--",
        linewidth=2,
    )

    ax.text(
        len(labels) - 0.2,
        THRESHOLD + 20,
        "1,000 SNP threshold",
        ha="right",
        va="bottom",
        fontsize=10,
        fontweight="bold",
    )


    # --------------------------------------------------------
    # VALUE LABELS
    # --------------------------------------------------------

    for bar, value, region in zip(
        bars,
        candidate_y,
        candidate_x,
    ):

        ax.annotate(
            f"{value:,}",
            xy=(
                bar.get_x()
                + bar.get_width() / 2,
                value,
            ),
            xytext=(0, 7),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=10,
            fontweight="bold",
        )


    # --------------------------------------------------------
    # TITLE
    # --------------------------------------------------------

    ax.set_title(
        "Chromosome 1 — Theta Regions Below 1,000 SNPs",
        fontsize=17,
        fontweight="bold",
        pad=15,
    )

    ax.grid(
        axis="y",
        linestyle="--",
        alpha=0.3,
    )

    ax.set_axisbelow(True)

    plt.tight_layout()

    candidate_graph_path = (
        OUTPUT_DIR
        / "chr1_theta_below_1000_zoom.png"
    )

    plt.savefig(
        candidate_graph_path,
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
    )

    plt.close()

else:
    candidate_graph_path = None


# ============================================================
# FINAL OUTPUT
# ============================================================

print()
print("=" * 75)
print("OUTPUT FILES")
print("=" * 75)

print(
    f"Overview graph:"
    f"\n  {overview_path}"
)

if candidate_graph_path:
    print(
        f"Candidate-only graph:"
        f"\n  {candidate_graph_path}"
    )

print(
    f"Candidate CSV:"
    f"\n  {STATS_DIR / 'chr1_regions_below_1000_snps.csv'}"
)

print()
print("Analysis completed successfully.")
print("=" * 75)