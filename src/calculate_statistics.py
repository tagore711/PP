import csv
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
INPUT = ROOT / "results" / "combined" / "all_chromosomes_region_counts.csv"
OUTPUT = ROOT / "results" / "statistics"
OUTPUT.mkdir(parents=True, exist_ok=True)

with INPUT.open(newline="", encoding="utf-8-sig") as f:
    regions = list(csv.DictReader(f))

for r in regions:
    for key in ("snp_count", "chromosome_length_bp", "theta_size_bp",
                "region_start_bp", "region_end_bp"):
        r[key] = int(r[key])

# Aggregate regional counts into chromosome totals
chromosomes = {}
for r in regions:
    c = r["chromosome"]
    chromosomes.setdefault(c, {
        "chromosome": c,
        "chromosome_length_bp": r["chromosome_length_bp"],
        "snp_count": 0,
        "number_of_regions": 0,
    })
    chromosomes[c]["snp_count"] += r["snp_count"]
    chromosomes[c]["number_of_regions"] += 1

chrom_list = sorted(
    chromosomes.values(),
    key=lambda x: int(x["chromosome"].replace("chr", ""))
)

total_snps = sum(c["snp_count"] for c in chrom_list)
counts = [c["snp_count"] for c in chrom_list]

# Regional statistics
region_counts = [r["snp_count"] for r in regions]

with (OUTPUT / "chromosome_statistics.csv").open(
    "w", newline="", encoding="utf-8"
) as f:
    fields = list(chrom_list[0].keys()) + [
        "snp_density_per_million_bp", "percentage_of_total_snps"
    ]
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    for c in chrom_list:
        c["snp_density_per_million_bp"] = round(
            c["snp_count"] / c["chromosome_length_bp"] * 1_000_000, 4
        )
        c["percentage_of_total_snps"] = round(
            c["snp_count"] / total_snps * 100, 4
        )
        writer.writerow(c)

overall = [
    ("chromosomes_processed", len(chrom_list)),
    ("theta_size_bp", 1_000_000),
    ("total_regions", len(regions)),
    ("total_snp_records", total_snps),
    ("mean_snps_per_chromosome", round(statistics.mean(counts), 2)),
    ("median_snps_per_chromosome", statistics.median(counts)),
    ("population_stddev_per_chromosome", round(statistics.pstdev(counts), 2)),
    ("minimum_snp_chromosome", min(chrom_list, key=lambda x: x["snp_count"])["chromosome"]),
    ("minimum_snp_count", min(counts)),
    ("maximum_snp_chromosome", max(chrom_list, key=lambda x: x["snp_count"])["chromosome"]),
    ("maximum_snp_count", max(counts)),
    ("mean_snps_per_region", round(statistics.mean(region_counts), 2)),
    ("median_snps_per_region", statistics.median(region_counts)),
    ("maximum_region_snp_count", max(region_counts)),
    ("minimum_region_snp_count", min(region_counts)),
]

with (OUTPUT / "overall_statistics.csv").open(
    "w", newline="", encoding="utf-8"
) as f:
    writer = csv.writer(f)
    writer.writerow(["metric", "value"])
    writer.writerows(overall)

print("\n=== OVERALL STATISTICS ===")
for k, v in overall:
    print(f"{k}: {v}")

print("\n=== TOP 5 CHROMOSOMES ===")
for c in sorted(chrom_list, key=lambda x: x["snp_count"], reverse=True)[:5]:
    print(f'{c["chromosome"]}: {c["snp_count"]:,}')

print("\nStatistics saved to:", OUTPUT)
