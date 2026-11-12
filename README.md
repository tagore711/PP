# Nucleotide-Transformer: Genomic θ-Region SNP Analysis

This repository contains the analysis pipeline, statistics, and visualization tools for studying single nucleotide polymorphism (SNP) density across 1 Mbp ($\theta$) genomic windows in human chromosomes (GRCh38).

## Overview

The pipeline analyzes whole-genome variant call format (VCF) data from the 1000 Genomes Project (phased biallelic SNVs and indels across chromosomes 1–22) to:
- Partition chromosomes into consecutive 1 Mbp ($\theta$) regions.
- Quantify SNP densities and identify regions with abnormally low variant frequencies ($\text{SNP count} < 10,000$).
- Generate per-chromosome and whole-genome summary statistics.
- Produce visualizations including heatmaps, chromosome-level distribution charts, and filtered region comparisons.

## Project Structure

```text
├── data/
│   └── raw/
│       ├── reference/
│       │   └── hg38.chrom.sizes      # GRCh38 chromosome size definitions
│       └── vcf/
│           └── README.md             # Instructions for 1000 Genomes VCF data
├── src/
│   ├── calculate_statistics.py       # Computes chromosome and overall statistics
│   ├── chr1_theta_filter_graph.py    # Generates filtered region plots for chr1
│   ├── combine_results.py            # Aggregates per-chromosome results
│   ├── count_snps.py                 # Core SNP counter per 1 Mbp theta window
│   ├── plot_snp_regions.py           # Plots SNP counts per genomic region
│   ├── plot_theta_below_10000.py     # Visualizations for regions with <10k SNPs
│   └── theta_split_heatmaps.py       # Heatmap generators for theta regions
├── results/
│   ├── combined/                     # Merged CSV summaries across all chromosomes
│   ├── per_chromosome/               # Per-chromosome counts and summaries (chr1-chr22)
│   ├── statistics/                   # Detailed statistical analyses
│   └── visualizations/               # Output plots and heatmaps
├── plot.py                           # Plotting utilities for low-SNP theta regions
├── theta_less.py                     # Identification & plotting of theta regions < 10,000 SNPs
└── README.md
```

## Requirements

- Python 3.9+
- `pandas`
- `numpy`
- `matplotlib`
- `seaborn`

To install dependencies:
```bash
pip install pandas numpy matplotlib seaborn
```

## Usage

1. **Download VCF Data**:
   Place the 1000 Genomes GRCh38 phased VCF files (`ALL.chr*.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz`) into `data/raw/vcf/`.

2. **Count SNPs**:
   ```bash
   python src/count_snps.py
   ```

3. **Combine Results & Generate Statistics**:
   ```bash
   python src/combine_results.py
   python src/calculate_statistics.py
   ```

4. **Generate Plots & Visualizations**:
   ```bash
   python plot.py
   python theta_less.py
   python src/theta_split_heatmaps.py
   ```
