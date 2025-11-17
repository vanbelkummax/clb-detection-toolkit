#!/usr/bin/env bash

###############################################################################
# Pipeline 1 Configuration
# Copy this to config.sh and customize for your environment
###############################################################################

# Executables (modify paths as needed)
BBMAP="/path/to/bbmap/bbmap.sh"
DIAMOND="/path/to/diamond"

# Reference databases
CLB_FASTA="/path/to/references/clb_complete_island.fasta"
CLB_PROTEIN_DB="/path/to/references/clb_complete_island.dmnd"

# Output directory
OUTDIR="./pipeline1_results"

# Parameters
THREADS=8
MIN_IDENTITY=0.90  # 90% nucleotide identity for BBMap
