# CLB Island Detection Toolkit

Comprehensive bioinformatics pipelines for detecting and characterizing colibactin (CLB) biosynthesis genes in metagenomic samples.

## Overview

This toolkit provides three complementary approaches for CLB island detection, each optimized for different research scenarios:

- **Pipeline 1 (bbmap)**: Read-based detection with high sensitivity for prevalence screening
- **Pipeline 2 (DIAMOND)**: Protein-level detection with host depletion for cleaner signal
- **Pipeline 3 (assembly)**: Structure validation and gene recovery from assembled genomes

## Features

- **Multi-method validation**: Cross-validate results across read-based and assembly-based approaches
- **Quantitative output**: RPM normalization, breadth/depth metrics, gene counts
- **Scalable**: Parallel processing for hundreds of samples
- **Flexible databases**: Support for custom CLB reference databases
- **Quality control**: Automated read QC, assembly quality checks, contamination filtering

## Quick Start

### Installation

```bash
# Clone repository
git clone https://github.com/YOUR_USERNAME/clb-detection-toolkit.git
cd clb-detection-toolkit

# Install dependencies (conda/mamba required)
bash utils/install_dependencies.sh

# Validate installation
bash utils/validate_installation.sh
```

### Running Pipelines

**Pipeline 1 - Read-based screening (recommended for prevalence studies):**
```bash
cd pipeline1_bbmap
cp config.example.sh config.sh
# Edit config.sh with your paths
bash run_pipeline1.sh sample_list.txt
```

**Pipeline 2 - DIAMOND protein search (recommended for low-biomass samples):**
```bash
cd pipeline2_diamond
cp config.example.sh config.sh
bash run_pipeline2.sh sample_list.txt
```

**Pipeline 3 - Assembly-based validation (recommended for complete island characterization):**
```bash
cd pipeline3_assembly
cp config.example.sh config.sh
bash run_pipeline3.sh assembly_list.txt
```

## Pipeline Comparison

| Feature | Pipeline 1 (bbmap) | Pipeline 2 (DIAMOND) | Pipeline 3 (assembly) |
|---------|-------------------|---------------------|----------------------|
| **Input** | Raw reads | Raw reads | Assembled contigs |
| **Sensitivity** | Very high | High | Moderate |
| **Specificity** | Moderate | High | Very high |
| **Host depletion** | No | Yes | N/A |
| **Gene structure** | Partial | Partial | Complete |
| **Runtime** | Fast (~5 min/sample) | Medium (~10 min/sample) | Slow (~20 min/sample) |
| **Best for** | Prevalence screening | Low-biomass samples | Taxonomic assignment |

## Output Files

All pipelines generate standardized output:

```
results/
├── stats/
│   ├── sample1_stats.json       # Per-sample metrics
│   └── summary_table.tsv         # Combined results table
├── hits/
│   ├── sample1_hits.tsv          # DIAMOND alignment results
│   └── sample1_genes_detected.txt
└── logs/
    └── sample1.log
```

## Requirements

- **Software**:
  - BBMap ≥38.0
  - Bowtie2 ≥2.4.0
  - DIAMOND ≥2.0.0
  - Prodigal ≥2.6.3
  - metaSPAdes ≥3.15.0 (for assembly)

- **Resources**:
  - RAM: 16 GB minimum, 32 GB recommended
  - Storage: ~500 MB per sample (fastq + intermediate files)
  - Threads: 8+ cores recommended

## Database Setup

Download pre-built CLB reference databases:

```bash
# Complete CLB island (19 genes, clbA-R)
wget https://example.com/clb_complete_island.fasta

# Build DIAMOND database
diamond makedb --in clb_complete_island.fasta --db clb_complete_island.dmnd
```

Or build custom database from GenBank accessions (see `docs/database_construction.md`).

## Citation

If you use this toolkit in your research, please consider listing me as a co-author and reach out to me
max.r.van.belkum@vanderbilt.edu
```

## License

MIT License - see LICENSE file for details.

## Contributing

Contributions welcome! Please see CONTRIBUTING.md for guidelines.


Developed for metagenomic analysis of colibactin-producing bacteria in infant gut microbiomes.

