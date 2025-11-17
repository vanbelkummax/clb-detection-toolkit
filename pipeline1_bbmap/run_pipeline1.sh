#!/usr/bin/env bash

###############################################################################
# Pipeline 1: Read-Based CLB Detection with BBMap
#
# High-sensitivity detection of CLB biosynthesis genes directly from
# metagenomic reads using BBMap alignment and DIAMOND BLASTX.
#
# Input:  Paired-end FASTQ files
# Output: Gene-level detection, RPM normalization, breadth/depth metrics
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
mkdir -p "$OUTDIR"/{bbmap_aligned,diamond_hits,stats,logs}

# Log configuration
echo "========================================" | tee "$OUTDIR/pipeline1.log"
echo "Pipeline 1: BBMap + DIAMOND BLASTX" | tee -a "$OUTDIR/pipeline1.log"
echo "========================================" | tee -a "$OUTDIR/pipeline1.log"
echo "Started: $(date)" | tee -a "$OUTDIR/pipeline1.log"
echo "Configuration:" | tee -a "$OUTDIR/pipeline1.log"
echo "  BBMap executable: $BBMAP" | tee -a "$OUTDIR/pipeline1.log"
echo "  DIAMOND executable: $DIAMOND" | tee -a "$OUTDIR/pipeline1.log"
echo "  CLB nucleotide ref: $CLB_FASTA" | tee -a "$OUTDIR/pipeline1.log"
echo "  CLB protein DB: $CLB_PROTEIN_DB" | tee -a "$OUTDIR/pipeline1.log"
echo "  Threads: $THREADS" | tee -a "$OUTDIR/pipeline1.log"
echo "" | tee -a "$OUTDIR/pipeline1.log"

# Process each sample
SAMPLE_COUNT=0
while IFS=$'\t' read -r sample_id read1 read2; do
    # Skip header
    [[ "$sample_id" == "SampleID" ]] && continue

    SAMPLE_COUNT=$((SAMPLE_COUNT + 1))
    echo ">>> Sample $SAMPLE_COUNT: $sample_id <<<" | tee -a "$OUTDIR/pipeline1.log"

    # Output files
    BBMAP_SAM="$OUTDIR/bbmap_aligned/${sample_id}_clb.sam"
    BBMAP_LOG="$OUTDIR/logs/${sample_id}_bbmap.log"
    UNALIGNED_FQ="$OUTDIR/bbmap_aligned/${sample_id}_unaligned.fastq"
    ALIGNED_FQ="$OUTDIR/bbmap_aligned/${sample_id}_aligned.fastq"
    DIAMOND_OUT="$OUTDIR/diamond_hits/${sample_id}_clb_hits.tsv"
    STATS_JSON="$OUTDIR/stats/${sample_id}_stats.json"

    # Step 1: BBMap alignment to CLB nucleotide reference
    echo "  [1/3] BBMap alignment..." | tee -a "$OUTDIR/pipeline1.log"
    $BBMAP \
        in1="$read1" \
        in2="$read2" \
        ref="$CLB_FASTA" \
        out="$BBMAP_SAM" \
        outu="$UNALIGNED_FQ" \
        outm="$ALIGNED_FQ" \
        minid="$MIN_IDENTITY" \
        threads="$THREADS" \
        nodisk \
        covstats="$OUTDIR/stats/${sample_id}_coverage.txt" \
        &> "$BBMAP_LOG"

    # Extract coverage statistics
    TOTAL_READS=$(grep "reads:" "$BBMAP_LOG" | head -1 | awk '{print $2}')
    MAPPED_READS=$(grep "mapped:" "$BBMAP_LOG" | head -1 | awk '{print $2}')

    echo "    Total reads: $TOTAL_READS" | tee -a "$OUTDIR/pipeline1.log"
    echo "    Mapped reads: $MAPPED_READS" | tee -a "$OUTDIR/pipeline1.log"

    # Step 2: DIAMOND BLASTX on aligned reads
    echo "  [2/3] DIAMOND BLASTX..." | tee -a "$OUTDIR/pipeline1.log"
    $DIAMOND blastx \
        --query "$ALIGNED_FQ" \
        --db "$CLB_PROTEIN_DB" \
        --out "$DIAMOND_OUT" \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
        --id "$MIN_IDENTITY" \
        --evalue 1e-10 \
        --threads "$THREADS" \
        --max-target-seqs 1 \
        2>&1 | tee -a "$OUTDIR/pipeline1.log"

    # Step 3: Calculate gene-level statistics
    echo "  [3/3] Calculating statistics..." | tee -a "$OUTDIR/pipeline1.log"

    # Count unique CLB genes detected
    if [ -f "$DIAMOND_OUT" ] && [ -s "$DIAMOND_OUT" ]; then
        GENES_DETECTED=$(grep -iE "clb[A-R]" "$DIAMOND_OUT" 2>/dev/null | cut -f2 | sort -u | wc -l || echo "0")
        TOTAL_HITS=$(wc -l < "$DIAMOND_OUT")
    else
        GENES_DETECTED=0
        TOTAL_HITS=0
    fi

    # Parse coverage for breadth/depth
    python3 << EOF > "$STATS_JSON"
import json
import sys

# Read coverage stats
coverage_data = {}
try:
    with open("$OUTDIR/stats/${sample_id}_coverage.txt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 8:
                gene = parts[0]
                avg_cov = float(parts[1])
                pct_covered = float(parts[5])
                coverage_data[gene] = {
                    "avg_depth": avg_cov,
                    "breadth": pct_covered / 100.0
                }
except:
    pass

# Calculate RPM
total_reads = int("$TOTAL_READS")
mapped_reads = int("$MAPPED_READS")
rpm = (mapped_reads / total_reads * 1e6) if total_reads > 0 else 0

# Output JSON
stats = {
    "sample_id": "$sample_id",
    "total_reads": total_reads,
    "mapped_reads": mapped_reads,
    "genes_detected": int("$GENES_DETECTED"),
    "total_hits": int("$TOTAL_HITS"),
    "clb_rpm": rpm,
    "coverage_data": coverage_data
}

print(json.dumps(stats, indent=2))
EOF

    echo "    Genes detected: $GENES_DETECTED" | tee -a "$OUTDIR/pipeline1.log"
    echo "    Total DIAMOND hits: $TOTAL_HITS" | tee -a "$OUTDIR/pipeline1.log"
    echo "" | tee -a "$OUTDIR/pipeline1.log"

done < "$SAMPLE_LIST"

# Generate summary table
echo "Generating summary table..." | tee -a "$OUTDIR/pipeline1.log"
python3 << 'EOF' > "$OUTDIR/stats/summary_table.tsv"
import json
import glob
import sys

# Header
print("SampleID\tTotalReads\tMappedReads\tGenes\tHits\tCLB_RPM\tMeanBreadth\tMeanDepth")

# Process all JSON files
for json_file in sorted(glob.glob("$OUTDIR/stats/*_stats.json")):
    with open(json_file) as f:
        data = json.load(f)

    # Calculate mean breadth/depth
    if data["coverage_data"]:
        breadths = [v["breadth"] for v in data["coverage_data"].values()]
        depths = [v["avg_depth"] for v in data["coverage_data"].values()]
        mean_breadth = sum(breadths) / len(breadths)
        mean_depth = sum(depths) / len(depths)
    else:
        mean_breadth = 0
        mean_depth = 0

    print(f"{data['sample_id']}\t{data['total_reads']}\t{data['mapped_reads']}\t"
          f"{data['genes_detected']}\t{data['total_hits']}\t{data['clb_rpm']:.2f}\t"
          f"{mean_breadth:.4f}\t{mean_depth:.2f}")
EOF

echo "" | tee -a "$OUTDIR/pipeline1.log"
echo "========================================" | tee -a "$OUTDIR/pipeline1.log"
echo "Pipeline 1 Complete!" | tee -a "$OUTDIR/pipeline1.log"
echo "Completed: $(date)" | tee -a "$OUTDIR/pipeline1.log"
echo "========================================" | tee -a "$OUTDIR/pipeline1.log"
echo "" | tee -a "$OUTDIR/pipeline1.log"
echo "Output files:" | tee -a "$OUTDIR/pipeline1.log"
echo "  Summary: $OUTDIR/stats/summary_table.tsv" | tee -a "$OUTDIR/pipeline1.log"
echo "  Per-sample stats: $OUTDIR/stats/" | tee -a "$OUTDIR/pipeline1.log"
echo "  DIAMOND hits: $OUTDIR/diamond_hits/" | tee -a "$OUTDIR/pipeline1.log"
echo "  Logs: $OUTDIR/logs/" | tee -a "$OUTDIR/pipeline1.log"
