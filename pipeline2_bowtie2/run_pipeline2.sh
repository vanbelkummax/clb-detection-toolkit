#!/usr/bin/env bash

###############################################################################
# Pipeline 2: Bowtie2 Host Depletion + DIAMOND Protein Detection
#
# Removes host contamination before CLB detection, ideal for low-biomass
# samples or host-associated microbiomes.
#
# Input:  Paired-end FASTQ files
# Output: Host-depleted reads, CLB gene detection with protein-level validation
###############################################################################

set -euo pipefail

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

# Usage
if [ $# -lt 1 ]; then
    echo "Usage: $0 <sample_list.txt>"
    echo ""
    echo "Sample list format (TSV):"
    echo "SampleID  Read1_Path  Read2_Path"
    exit 1
fi

SAMPLE_LIST=$1

# Create output directories
mkdir -p "$OUTDIR"/{qc_reads,host_depleted,diamond_hits,stats,logs}

# Log configuration
echo "========================================" | tee "$OUTDIR/pipeline2.log"
echo "Pipeline 2: Bowtie2 + DIAMOND" | tee -a "$OUTDIR/pipeline2.log"
echo "========================================" | tee -a "$OUTDIR/pipeline2.log"
echo "Started: $(date)" | tee -a "$OUTDIR/pipeline2.log"
echo "Configuration:" | tee -a "$OUTDIR/pipeline2.log"
echo "  bbduk: $BBDUK" | tee -a "$OUTDIR/pipeline2.log"
echo "  Bowtie2: $BOWTIE2" | tee -a "$OUTDIR/pipeline2.log"
echo "  DIAMOND: $DIAMOND" | tee -a "$OUTDIR/pipeline2.log"
echo "  Host reference: $HOST_REFERENCE" | tee -a "$OUTDIR/pipeline2.log"
echo "  CLB protein DB: $CLB_PROTEIN_DB" | tee -a "$OUTDIR/pipeline2.log"
echo "  Threads: $THREADS" | tee -a "$OUTDIR/pipeline2.log"
echo "" | tee -a "$OUTDIR/pipeline2.log"

# Process each sample
SAMPLE_COUNT=0
while IFS=$'\t' read -r sample_id read1 read2; do
    # Skip header
    [[ "$sample_id" == "SampleID" ]] && continue

    SAMPLE_COUNT=$((SAMPLE_COUNT + 1))
    echo ">>> Sample $SAMPLE_COUNT: $sample_id <<<" | tee -a "$OUTDIR/pipeline2.log"

    # Output files
    QC_R1="$OUTDIR/qc_reads/${sample_id}_R1_qc.fastq.gz"
    QC_R2="$OUTDIR/qc_reads/${sample_id}_R2_qc.fastq.gz"
    DEPLETED_R1="$OUTDIR/host_depleted/${sample_id}_R1_depleted.fastq.gz"
    DEPLETED_R2="$OUTDIR/host_depleted/${sample_id}_R2_depleted.fastq.gz"
    DIAMOND_OUT="$OUTDIR/diamond_hits/${sample_id}_clb_hits.tsv"
    STATS_JSON="$OUTDIR/stats/${sample_id}_stats.json"

    # Step 1: Quality trimming with bbduk
    echo "  [1/4] Quality trimming..." | tee -a "$OUTDIR/pipeline2.log"
    $BBDUK \
        in1="$read1" \
        in2="$read2" \
        out1="$QC_R1" \
        out2="$QC_R2" \
        ref="$ADAPTERS" \
        ktrim=r k=23 mink=11 hdist=1 \
        qtrim=rl trimq=20 \
        minlength=50 \
        threads="$THREADS" \
        &> "$OUTDIR/logs/${sample_id}_bbduk.log"

    READS_AFTER_QC=$(zcat "$QC_R1" | wc -l | awk '{print $1/4}')
    echo "    Reads after QC: $READS_AFTER_QC" | tee -a "$OUTDIR/pipeline2.log"

    # Step 2: Host depletion with Bowtie2
    echo "  [2/4] Host depletion..." | tee -a "$OUTDIR/pipeline2.log"
    $BOWTIE2 \
        -x "$HOST_REFERENCE" \
        -1 "$QC_R1" \
        -2 "$QC_R2" \
        --un-conc-gz "$OUTDIR/host_depleted/${sample_id}_R%_depleted.fastq.gz" \
        --threads "$THREADS" \
        --no-unal \
        2> "$OUTDIR/logs/${sample_id}_bowtie2.log" \
        > /dev/null

    READS_AFTER_DEPLETION=$(zcat "$DEPLETED_R1" | wc -l | awk '{print $1/4}')
    HOST_READS=$((READS_AFTER_QC - READS_AFTER_DEPLETION))
    echo "    Reads after host depletion: $READS_AFTER_DEPLETION" | tee -a "$OUTDIR/pipeline2.log"
    echo "    Host reads removed: $HOST_READS" | tee -a "$OUTDIR/pipeline2.log"

    # Step 3: Interleave reads for DIAMOND
    echo "  [3/4] DIAMOND BLASTX..." | tee -a "$OUTDIR/pipeline2.log"
    INTERLEAVED="$OUTDIR/host_depleted/${sample_id}_interleaved.fastq.gz"
    paste <(zcat "$DEPLETED_R1") <(zcat "$DEPLETED_R2") | \
        awk '{print $1"\n"$2"\n"$3"\n"$4}' | \
        gzip > "$INTERLEAVED"

    $DIAMOND blastx \
        --query "$INTERLEAVED" \
        --db "$CLB_PROTEIN_DB" \
        --out "$DIAMOND_OUT" \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
        --id "$MIN_IDENTITY" \
        --evalue 1e-10 \
        --threads "$THREADS" \
        --max-target-seqs 1 \
        2>&1 | tee -a "$OUTDIR/logs/${sample_id}_diamond.log"

    # Step 4: Calculate statistics
    echo "  [4/4] Calculating statistics..." | tee -a "$OUTDIR/pipeline2.log"

    if [ -f "$DIAMOND_OUT" ] && [ -s "$DIAMOND_OUT" ]; then
        GENES_DETECTED=$(grep -iE "clb[A-R]" "$DIAMOND_OUT" 2>/dev/null | cut -f2 | sort -u | wc -l || echo "0")
        TOTAL_HITS=$(wc -l < "$DIAMOND_OUT")
    else
        GENES_DETECTED=0
        TOTAL_HITS=0
    fi

    # Calculate RPM based on post-depletion reads
    RPM=$(awk "BEGIN {print ($TOTAL_HITS / $READS_AFTER_DEPLETION * 1e6)}")

    # Generate JSON stats
    cat > "$STATS_JSON" << EOF
{
  "sample_id": "$sample_id",
  "reads_after_qc": $READS_AFTER_QC,
  "host_reads_removed": $HOST_READS,
  "reads_after_depletion": $READS_AFTER_DEPLETION,
  "genes_detected": $GENES_DETECTED,
  "total_hits": $TOTAL_HITS,
  "clb_rpm": $RPM
}
EOF

    echo "    Genes detected: $GENES_DETECTED" | tee -a "$OUTDIR/pipeline2.log"
    echo "    CLB RPM: $RPM" | tee -a "$OUTDIR/pipeline2.log"
    echo "" | tee -a "$OUTDIR/pipeline2.log"

    # Cleanup interleaved file
    rm -f "$INTERLEAVED"

done < "$SAMPLE_LIST"

# Generate summary table
echo "Generating summary table..." | tee -a "$OUTDIR/pipeline2.log"
{
    echo -e "SampleID\tReadsAfterQC\tHostReads\tReadsAfterDepletion\tGenes\tHits\tCLB_RPM"
    for json_file in "$OUTDIR/stats"/*_stats.json; do
        python3 -c "
import json
with open('$json_file') as f:
    d = json.load(f)
print(f\"{d['sample_id']}\t{d['reads_after_qc']}\t{d['host_reads_removed']}\t{d['reads_after_depletion']}\t{d['genes_detected']}\t{d['total_hits']}\t{d['clb_rpm']:.2f}\")
"
    done
} > "$OUTDIR/stats/summary_table.tsv"

echo "" | tee -a "$OUTDIR/pipeline2.log"
echo "========================================" | tee -a "$OUTDIR/pipeline2.log"
echo "Pipeline 2 Complete!" | tee -a "$OUTDIR/pipeline2.log"
echo "Completed: $(date)" | tee -a "$OUTDIR/pipeline2.log"
echo "========================================" | tee -a "$OUTDIR/pipeline2.log"
echo "" | tee -a "$OUTDIR/pipeline2.log"
echo "Output files:" | tee -a "$OUTDIR/pipeline2.log"
echo "  Summary: $OUTDIR/stats/summary_table.tsv" | tee -a "$OUTDIR/pipeline2.log"
echo "  Host-depleted reads: $OUTDIR/host_depleted/" | tee -a "$OUTDIR/pipeline2.log"
echo "  DIAMOND hits: $OUTDIR/diamond_hits/" | tee -a "$OUTDIR/pipeline2.log"
