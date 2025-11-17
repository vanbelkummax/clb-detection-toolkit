#!/usr/bin/env bash

###############################################################################
# Pipeline 3 Configuration
# Copy this to config.sh and customize for your environment
###############################################################################

# Executables
PRODIGAL="/path/to/prodigal"
DIAMOND="/path/to/diamond"

# Reference databases
CLB_PROTEIN_DB="/path/to/references/clb_complete_island.dmnd"

# Output directory
OUTDIR="./pipeline3_results"

# Parameters
THREADS=8
MIN_IDENTITY=0.90  # 90% amino acid identity for DIAMOND BLASTP
