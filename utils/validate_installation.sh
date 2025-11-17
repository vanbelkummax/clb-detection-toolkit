#!/usr/bin/env bash

###############################################################################
# Installation Validation Script
# Checks that all required tools are available
###############################################################################

set -euo pipefail

echo "========================================="
echo "Validating CLB Detection Toolkit Installation"
echo "========================================="
echo ""

ERRORS=0

# Function to check command
check_command() {
    local cmd=$1
    local name=$2

    if command -v "$cmd" &> /dev/null; then
        local version=$($cmd --version 2>&1 | head -1 || echo "unknown")
        echo "✓ $name: $version"
    else
        echo "✗ $name: NOT FOUND"
        ERRORS=$((ERRORS + 1))
    fi
}

# Check core tools
echo "Checking required tools..."
echo ""

check_command "bbmap.sh" "BBMap"
check_command "bowtie2" "Bowtie2"
check_command "diamond" "DIAMOND"
check_command "prodigal" "Prodigal"
check_command "metaspades.py" "metaSPAdes"
check_command "samtools" "Samtools"
check_command "python3" "Python3"

echo ""
echo "========================================="
if [ $ERRORS -eq 0 ]; then
    echo "✓ All required tools are installed!"
    echo "========================================="
    echo ""
    echo "Next steps:"
    echo "1. Download or build CLB reference databases"
    echo "2. Copy config.example.sh to config.sh in each pipeline directory"
    echo "3. Edit config.sh with your paths"
    echo "4. Run pipelines!"
    exit 0
else
    echo "✗ Missing $ERRORS required tool(s)"
    echo "========================================="
    echo ""
    echo "Please run: bash utils/install_dependencies.sh"
    exit 1
fi
