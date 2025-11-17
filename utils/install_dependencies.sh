#!/usr/bin/env bash

###############################################################################
# Dependency Installation Script
# Installs all required tools using conda/mamba
###############################################################################

set -euo pipefail

echo "========================================="
echo "CLB Detection Toolkit - Dependency Setup"
echo "========================================="
echo ""

# Check for conda/mamba
if command -v mamba &> /dev/null; then
    CONDA_CMD="mamba"
    echo "✓ Using mamba for installation"
elif command -v conda &> /dev/null; then
    CONDA_CMD="conda"
    echo "✓ Using conda for installation"
else
    echo "ERROR: Neither conda nor mamba found!"
    echo "Please install miniconda or miniforge first:"
    echo "  https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

echo ""
echo "Creating conda environment: clb_detection"
echo ""

# Create environment with all dependencies
$CONDA_CMD create -n clb_detection -y \
    python=3.9 \
    bbmap \
    bowtie2 \
    diamond \
    prodigal \
    spades \
    samtools \
    seqtk \
    -c bioconda -c conda-forge

echo ""
echo "========================================="
echo "Installation Complete!"
echo "========================================="
echo ""
echo "To activate the environment:"
echo "  conda activate clb_detection"
echo ""
echo "To verify installation:"
echo "  bash utils/validate_installation.sh"
echo ""
