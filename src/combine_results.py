import csv
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
INPUT_DIR = ROOT / "results" / "per_chromosome"
OUTPUT_DIR = ROOT / "results" / "combined"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

all_regions = []
chromosome_summary = []

# Validate and collect results for chromosomes 1–22
for chrom_num in range(1, 23):
    chrom = f"chr{chrom_num}"
    folder = INPUT_DIR / chrom

    region_file = folder / "region_counts.csv"
    summary_file = folder / "processing_summary.csv"

    if not region_file.exists() or not summary_file.exists():
        raise FileNotFoundError(f"Missing output file for {chrom}")

    with region_file.open("r", newline="", encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))

    if not rows:
        raise ValueError(f"No region data found in {region_file}")

    for row in rows:
        row["chromosome"] = chrom
        all_regions.append(row)

    with summary_file.open("r", newline="", encoding="utf-8-sig") as f:
        summary_rows = list(csv.DictReader(f))

    if not summary_rows:
        raise ValueError(f"No summary data found in {summary_file}")

    summary = summary_rows[0]
    summary["chromosome"] = chrom
    chromosome_summary.append(summary)

# Write combined region counts
region_path = OUTPUT_DIR / "all_chromosomes_region_counts.csv"
with region_path.open("w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=["chromosome"] + [
        k for k in all_regions[0] if k != "chromosome"
    ])
    writer.writeheader()
    writer.writerows(all_regions)

# Write chromosome-level processing summaries
summary_path = OUTPUT_DIR / "chromosome_summary.csv"
with summary_path.open("w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=["chromosome"] + [
            k for k in chromosome_summary[0] if k != "chromosome"
        ]
    )
    writer.writeheader()
    writer.writerows(chromosome_summary)

print("COMBINATION COMPLETED")
print(f"Chromosomes processed: {len(chromosome_summary)}")
print(f"Region rows combined: {len(all_regions)}")
print(f"Combined region counts: {region_path}")
print(f"Chromosome summary: {summary_path}")
