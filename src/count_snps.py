
import argparse
import csv
import gzip
import math
import re
from collections import defaultdict
from pathlib import Path


# ============================================================
# CONFIGURATION
# ============================================================

PROJECT_DIR = Path(__file__).resolve().parent.parent

VCF_DIR = PROJECT_DIR / "data" / "raw" / "vcf"
REFERENCE_DIR = PROJECT_DIR / "data" / "raw" / "reference"
RESULTS_DIR = PROJECT_DIR / "results"

CHROM_SIZES_FILE = REFERENCE_DIR / "hg38.chrom.sizes"

THETA_SIZE_BP = 1_000_000  # 1 million bp (1 Mbp)


# ============================================================
# LOAD GRCh38 CHROMOSOME LENGTHS
# ============================================================

def load_chromosome_sizes():
    """Load chromosome lengths from the GRCh38 reference file."""

    if not CHROM_SIZES_FILE.exists():
        raise FileNotFoundError(
            f"Chromosome sizes file not found: {CHROM_SIZES_FILE}\n"
            "Download hg38.chrom.sizes before running the script."
        )

    sizes = {}

    with open(CHROM_SIZES_FILE, "r", encoding="utf-8") as f:
        for line in f:
            fields = line.strip().split()

            if len(fields) < 2:
                continue

            chrom = fields[0].removeprefix("chr")
            length = int(fields[1])

            sizes[chrom] = length

    if not sizes:
        raise ValueError("No chromosome lengths were loaded.")

    return sizes


# ============================================================
# FILE HANDLING
# ============================================================

def open_vcf(path):
    """Open plain-text or gzip-compressed VCF files."""

    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")

    return open(path, "rt", encoding="utf-8")


def parse_contig_length(line):
    """Extract chromosome ID and length from a VCF contig header."""

    match = re.search(r"##contig=<(.+)>", line)

    if not match:
        return None, None

    attributes = dict(
        re.findall(r'(\w+)=("[^"]*"|[^,>]+)', match.group(1))
    )

    chrom = attributes.get("ID")
    length = attributes.get("length")

    if chrom:
        chrom = chrom.removeprefix("chr")

    if chrom and length and length.isdigit():
        return chrom, int(length)

    return chrom, None


# ============================================================
# SNP CLASSIFICATION
# ============================================================

def is_snp(ref, alt):
    """
    Return True for single-base SNP records.

    Accepts A/C/G/T REF and ALT alleles.
    Multiallelic records are accepted only when every ALT
    allele is a single-base substitution.
    """

    valid_bases = {"A", "C", "G", "T"}

    if ref not in valid_bases:
        return False

    alt_alleles = alt.split(",")

    return all(
        allele in valid_bases and len(allele) == 1
        for allele in alt_alleles
    )


# ============================================================
# PROCESS ONE VCF FILE FOR ONE CHROMOSOME
# ============================================================

def process_vcf(vcf_path, reference_sizes, target_chrom):
    print(f"\nProcessing: {vcf_path.name}")
    print(f"Target chromosome: chr{target_chrom}")

    # SNP counts indexed by chromosome and region index
    region_snp_counts = defaultdict(int)

    total_records = defaultdict(int)
    excluded_records = defaultdict(int)

    header_lengths = {}
    max_positions = defaultdict(int)

    # --------------------------------------------------------
    # Read VCF
    # --------------------------------------------------------

    with open_vcf(vcf_path) as vcf:

        for line in vcf:

            if line.startswith("##contig="):
                chrom, length = parse_contig_length(line)

                if chrom == target_chrom and length:
                    header_lengths[chrom] = length

                continue

            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 8:
                continue

            chrom = fields[0].removeprefix("chr")

            # Process only the selected chromosome.
            if chrom != target_chrom:
                continue

            pos = int(fields[1])
            ref = fields[3].upper()
            alt = fields[4].upper()

            total_records[chrom] += 1

            max_positions[chrom] = max(
                max_positions[chrom], pos
            )

            if pos < 1:
                continue

            # Region index is zero-based.
            # Position 1 belongs to region index 0.
            region_index = (pos - 1) // THETA_SIZE_BP

            if is_snp(ref, alt):
                region_snp_counts[(chrom, region_index)] += 1
            else:
                excluded_records[chrom] += 1

    # --------------------------------------------------------
    # Build region-level output
    # --------------------------------------------------------

    region_rows = []
    detail_rows = []

    all_chromosomes = (
        set(total_records)
        | set(header_lengths)
    )

    # Only the requested chromosome is included.
    all_chromosomes = {
        chrom for chrom in all_chromosomes
        if chrom == target_chrom
    }

    for chrom in sorted(
        all_chromosomes,
        key=lambda x: (
            0, int(x)
        ) if x.isdigit() else (1, x)
    ):

        # Prefer the official GRCh38 reference length.
        if chrom in reference_sizes:
            chromosome_length = reference_sizes[chrom]
            length_source = "GRCh38 reference"

        elif chrom in header_lengths:
            chromosome_length = header_lengths[chrom]
            length_source = "VCF header"

        else:
            chromosome_length = max_positions[chrom]
            length_source = "Maximum observed position (fallback)"

        if chromosome_length <= 0:
            continue

        # Number of theta-regions
        number_of_regions = math.ceil(
            chromosome_length / THETA_SIZE_BP
        )

        chromosome_snp_total = 0

        print(
            f"  Chromosome {chrom}: "
            f"{number_of_regions} theta-region(s)"
        )

        # ----------------------------------------------------
        # Generate one row per region
        # ----------------------------------------------------

        for region_index in range(number_of_regions):

            start_bp = (
                region_index * THETA_SIZE_BP
            ) + 1

            end_bp = min(
                (region_index + 1) * THETA_SIZE_BP,
                chromosome_length
            )

            snp_count = region_snp_counts.get(
                (chrom, region_index), 0
            )

            chromosome_snp_total += snp_count

            region_rows.append({
                "chromosome": chrom,
                "chromosome_length_bp": chromosome_length,
                "length_source": length_source,
                "theta_size_bp": THETA_SIZE_BP,
                "number_of_theta_regions": number_of_regions,
                "region_number": region_index + 1,
                "region_start_bp": start_bp,
                "region_end_bp": end_bp,
                "snp_count": snp_count,
            })

            print(
                f"    Region {region_index + 1}: "
                f"{start_bp:,} - {end_bp:,} bp | "
                f"SNPs: {snp_count:,}"
            )

        # ----------------------------------------------------
        # File processing details
        # ----------------------------------------------------

        detail_rows.append({
            "vcf_file": vcf_path.name,
            "chromosome": chrom,
            "chromosome_length_bp": chromosome_length,
            "length_source": length_source,
            "number_of_theta_regions": number_of_regions,
            "total_vcf_records": total_records[chrom],
            "snp_records": chromosome_snp_total,
            "excluded_records": excluded_records[chrom],
        })

    return region_rows, detail_rows


# ============================================================
# WRITE CSV
# ============================================================

def write_csv(output_path, rows):

    if not rows:
        print(f"No data to write: {output_path}")
        return

    with open(
        output_path,
        "w",
        newline="",
        encoding="utf-8"
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=list(rows[0].keys())
        )

        writer.writeheader()
        writer.writerows(rows)

    print(f"Saved: {output_path}")


# ============================================================
# MAIN PIPELINE
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description="Count SNP records for one chromosome."
    )

    parser.add_argument(
        "--chromosome",
        required=True,
        type=int,
        choices=range(1, 23),
        help="Autosome number from 1 to 22."
    )

    args = parser.parse_args()

    chrom = str(args.chromosome)
    chrom_label = f"chr{chrom}"

    # --------------------------------------------------------
    # Create chromosome-specific output directory
    # --------------------------------------------------------

    output_dir = RESULTS_DIR / "per_chromosome" / chrom_label
    output_dir.mkdir(parents=True, exist_ok=True)

    summary_file = output_dir / "region_counts.csv"
    details_file = output_dir / "processing_summary.csv"

    print("=" * 60)
    print(f"SNP COUNTING PIPELINE â€” {chrom_label}")
    print("=" * 60)

    print(f"Theta size: {THETA_SIZE_BP:,} bp")

    # --------------------------------------------------------
    # Load reference chromosome sizes
    # --------------------------------------------------------

    reference_sizes = load_chromosome_sizes()

    if chrom not in reference_sizes:
        raise ValueError(
            f"Chromosome {chrom} not found in reference sizes."
        )

    print(
        f"Reference chromosome length: "
        f"{reference_sizes[chrom]:,} bp"
    )

    # --------------------------------------------------------
    # Find only the selected chromosome's VCF file
    # --------------------------------------------------------

    vcf_files = (
        list(VCF_DIR.glob(f"*.chr{chrom}.*.vcf.gz"))
        + list(VCF_DIR.glob(f"*.chr{chrom}.*.vcf"))
    )

    if len(vcf_files) != 1:
        raise FileNotFoundError(
            f"Expected exactly one VCF for {chrom_label}, "
            f"but found {len(vcf_files)}: {vcf_files}"
        )

    vcf_path = vcf_files[0]

    print(f"Selected VCF: {vcf_path.name}")

    # --------------------------------------------------------
    # Process selected chromosome
    # --------------------------------------------------------

    region_rows, detail_rows = process_vcf(
        vcf_path,
        reference_sizes,
        chrom
    )

    if not detail_rows:
        raise ValueError(
            f"No records found for {chrom_label}. "
            "Check the VCF chromosome naming and file contents."
        )

    # --------------------------------------------------------
    # Save results
    # --------------------------------------------------------

    write_csv(summary_file, region_rows)
    write_csv(details_file, detail_rows)

    # --------------------------------------------------------
    # Final summary
    # --------------------------------------------------------

    total_snps = sum(
        row["snp_count"] for row in region_rows
    )

    print("\n" + "=" * 60)
    print(f"{chrom_label} PROCESSING COMPLETED")
    print("=" * 60)

    print(f"Region summary: {summary_file}")
    print(f"Processing details: {details_file}")
    print(f"Total region rows: {len(region_rows)}")
    print(f"Total SNP records: {total_snps:,}")


if __name__ == "__main__":
    main()
