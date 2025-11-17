#!/usr/bin/env bash

###############################################################################
# Pipeline 3: Assembly-Based CLB Detection
#
# Validates CLB island structure and enables gene recovery from assembled
# metagenomes using Prodigal gene prediction and DIAMOND BLASTP.
#
# Input:  Assembled contigs/scaffolds (FASTA)
# Output: Complete gene detection, scaffold-level annotations, gene sequences
###############################################################################

set -euo pipefail

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

# Usage
if [ $# -lt 1 ]; then
    echo "Usage: $0 <assembly_list.txt>"
    echo ""
    echo "Assembly list format (TSV):"
    echo "SampleID  Assembly_Path"
    exit 1
fi

ASSEMBLY_LIST=$1

# Create output directories
mkdir -p "$OUTDIR"/{prodigal_proteins,prodigal_genes,diamond_hits,extracted_sequences,stats,logs}

# Log configuration
echo "========================================" | tee "$OUTDIR/pipeline3.log"
echo "Pipeline 3: Assembly-Based Detection" | tee -a "$OUTDIR/pipeline3.log"
echo "========================================" | tee -a "$OUTDIR/pipeline3.log"
echo "Started: $(date)" | tee -a "$OUTDIR/pipeline3.log"
echo "Configuration:" | tee -a "$OUTDIR/pipeline3.log"
echo "  Prodigal: $PRODIGAL" | tee -a "$OUTDIR/pipeline3.log"
echo "  DIAMOND: $DIAMOND" | tee -a "$OUTDIR/pipeline3.log"
echo "  CLB protein DB: $CLB_PROTEIN_DB" | tee -a "$OUTDIR/pipeline3.log"
echo "  Threads: $THREADS" | tee -a "$OUTDIR/pipeline3.log"
echo "" | tee -a "$OUTDIR/pipeline3.log"

# Process each assembly
SAMPLE_COUNT=0
while IFS=$'\t' read -r sample_id assembly_path; do
    # Skip header
    [[ "$sample_id" == "SampleID" ]] && continue

    SAMPLE_COUNT=$((SAMPLE_COUNT + 1))
    echo ">>> Sample $SAMPLE_COUNT: $sample_id <<<" | tee -a "$OUTDIR/pipeline3.log"

    # Check if assembly exists
    if [ ! -f "$assembly_path" ]; then
        echo "  ERROR: Assembly not found: $assembly_path" | tee -a "$OUTDIR/pipeline3.log"
        echo "{\"sample_id\": \"$sample_id\", \"error\": \"Assembly not found\"}" > "$OUTDIR/stats/${sample_id}_stats.json"
        continue
    fi

    # Output files
    PROTEINS="$OUTDIR/prodigal_proteins/${sample_id}_proteins.faa"
    GENES="$OUTDIR/prodigal_genes/${sample_id}_genes.fna"
    DIAMOND_OUT="$OUTDIR/diamond_hits/${sample_id}_clb_hits.tsv"
    STATS_JSON="$OUTDIR/stats/${sample_id}_stats.json"

    # Step 1: Gene prediction with Prodigal
    echo "  [1/4] Running Prodigal..." | tee -a "$OUTDIR/pipeline3.log"
    $PRODIGAL \
        -i "$assembly_path" \
        -a "$PROTEINS" \
        -d "$GENES" \
        -p meta \
        &> "$OUTDIR/logs/${sample_id}_prodigal.log"

    NUM_PROTEINS=$(grep -c "^>" "$PROTEINS" 2>/dev/null || echo "0")
    echo "    Predicted proteins: $NUM_PROTEINS" | tee -a "$OUTDIR/pipeline3.log"

    if [ "$NUM_PROTEINS" -eq 0 ]; then
        echo "  WARNING: No proteins predicted!" | tee -a "$OUTDIR/pipeline3.log"
        echo "{\"sample_id\": \"$sample_id\", \"proteins_predicted\": 0, \"genes_detected\": 0}" > "$STATS_JSON"
        continue
    fi

    # Step 2: DIAMOND BLASTP protein vs protein
    echo "  [2/4] Running DIAMOND BLASTP..." | tee -a "$OUTDIR/pipeline3.log"
    $DIAMOND blastp \
        --query "$PROTEINS" \
        --db "$CLB_PROTEIN_DB" \
        --out "$DIAMOND_OUT" \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
        --id "$MIN_IDENTITY" \
        --evalue 1e-10 \
        --threads "$THREADS" \
        --max-target-seqs 1 \
        &> "$OUTDIR/logs/${sample_id}_diamond.log"

    # Step 3: Count CLB genes and scaffolds
    echo "  [3/4] Analyzing results..." | tee -a "$OUTDIR/pipeline3.log"

    if [ -f "$DIAMOND_OUT" ] && [ -s "$DIAMOND_OUT" ]; then
        # Count unique CLB genes (clbA-R)
        GENES_DETECTED=$(grep -iE "clb[A-R]" "$DIAMOND_OUT" 2>/dev/null | cut -f2 | sort -u | wc -l || echo "0")

        # Count scaffolds with CLB hits
        SCAFFOLDS_WITH_CLB=$(grep -iE "clb[A-R]" "$DIAMOND_OUT" 2>/dev/null | cut -f1 | cut -d'_' -f1-4 | sort -u | wc -l || echo "0")

        # Total CLB hits
        TOTAL_HITS=$(grep -icE "clb[A-R]" "$DIAMOND_OUT" 2>/dev/null || echo "0")

        # List detected genes
        DETECTED_GENES=$(grep -iE "clb[A-R]" "$DIAMOND_OUT" 2>/dev/null | cut -f2 | sort -u | tr '\n' ',' | sed 's/,$//' || echo "")
    else
        GENES_DETECTED=0
        SCAFFOLDS_WITH_CLB=0
        TOTAL_HITS=0
        DETECTED_GENES=""
    fi

    echo "    CLB genes detected: $GENES_DETECTED" | tee -a "$OUTDIR/pipeline3.log"
    echo "    CLB-containing scaffolds: $SCAFFOLDS_WITH_CLB" | tee -a "$OUTDIR/pipeline3.log"
    echo "    Total hits: $TOTAL_HITS" | tee -a "$OUTDIR/pipeline3.log"

    # Step 4: Extract CLB gene sequences
    if [ "$GENES_DETECTED" -gt 0 ]; then
        echo "  [4/4] Extracting CLB sequences..." | tee -a "$OUTDIR/pipeline3.log"

        EXTRACTED_PROTEINS="$OUTDIR/extracted_sequences/${sample_id}_clb_proteins.faa"
        EXTRACTED_GENES="$OUTDIR/extracted_sequences/${sample_id}_clb_genes.fna"

        # Extract protein sequences
        grep -iE "clb[A-R]" "$DIAMOND_OUT" | cut -f1 | sort -u | while read gene_id; do
            grep -A1 "^>$gene_id" "$PROTEINS"
        done > "$EXTRACTED_PROTEINS"

        # Extract nucleotide sequences
        grep -iE "clb[A-R]" "$DIAMOND_OUT" | cut -f1 | sort -u | while read gene_id; do
            grep -A1 "^>$gene_id" "$GENES"
        done > "$EXTRACTED_GENES"

        echo "    Extracted sequences saved" | tee -a "$OUTDIR/pipeline3.log"
    fi

    # Generate JSON stats
    cat > "$STATS_JSON" << EOF
{
  "sample_id": "$sample_id",
  "proteins_predicted": $NUM_PROTEINS,
  "genes_detected": $GENES_DETECTED,
  "scaffolds_with_clb": $SCAFFOLDS_WITH_CLB,
  "total_hits": $TOTAL_HITS,
  "detected_genes": "$DETECTED_GENES",
  "island_completeness": "$(awk "BEGIN {printf \"%.1f\", $GENES_DETECTED/19*100}")%"
}
EOF

    echo "" | tee -a "$OUTDIR/pipeline3.log"

done < "$ASSEMBLY_LIST"

# Generate summary table
echo "Generating summary table..." | tee -a "$OUTDIR/pipeline3.log"
{
    echo -e "SampleID\tProteins\tGenes\tScaffolds\tHits\tCompleteness\tDetectedGenes"
    for json_file in "$OUTDIR/stats"/*_stats.json; do
        python3 -c "
import json
import sys
try:
    with open('$json_file') as f:
        d = json.load(f)
    if 'error' in d:
        print(f\"{d['sample_id']}\tNA\tNA\tNA\tNA\tNA\t{d.get('error', 'Unknown error')}\")
    else:
        print(f\"{d['sample_id']}\t{d.get('proteins_predicted', 0)}\t{d.get('genes_detected', 0)}\t{d.get('scaffolds_with_clb', 0)}\t{d.get('total_hits', 0)}\t{d.get('island_completeness', '0%')}\t{d.get('detected_genes', '')}\")
except Exception as e:
    print(f\"Error processing $json_file: {e}\", file=sys.stderr)
"
    done
} > "$OUTDIR/stats/summary_table.tsv"

echo "" | tee -a "$OUTDIR/pipeline3.log"
echo "========================================" | tee -a "$OUTDIR/pipeline3.log"
echo "Pipeline 3 Complete!" | tee -a "$OUTDIR/pipeline3.log"
echo "Completed: $(date)" | tee -a "$OUTDIR/pipeline3.log"
echo "========================================" | tee -a "$OUTDIR/pipeline3.log"
echo "" | tee -a "$OUTDIR/pipeline3.log"
echo "Output files:" | tee -a "$OUTDIR/pipeline3.log"
echo "  Summary: $OUTDIR/stats/summary_table.tsv" | tee -a "$OUTDIR/pipeline3.log"
echo "  Predicted proteins: $OUTDIR/prodigal_proteins/" | tee -a "$OUTDIR/pipeline3.log"
echo "  Extracted CLB sequences: $OUTDIR/extracted_sequences/" | tee -a "$OUTDIR/pipeline3.log"
echo "  DIAMOND hits: $OUTDIR/diamond_hits/" | tee -a "$OUTDIR/pipeline3.log"
