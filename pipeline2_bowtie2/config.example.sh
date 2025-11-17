#!/usr/bin/env bash

###############################################################################
# Pipeline 2 Configuration
# Copy this to config.sh and customize for your environment
###############################################################################

# Executables
BBDUK="/path/to/bbmap/bbduk.sh"
BOWTIE2="/path/to/bowtie2"
DIAMOND="/path/to/diamond"

# Reference databases
HOST_REFERENCE="/path/to/host_reference/grch38"  # Bowtie2 index basename
ADAPTERS="/path/to/bbmap/resources/adapters.fa"
CLB_PROTEIN_DB="/path/to/references/clb_complete_island.dmnd"

# Output directory
OUTDIR="./pipeline2_results"

# Parameters
THREADS=8
MIN_IDENTITY=0.90  # 90% amino acid identity for DIAMOND
