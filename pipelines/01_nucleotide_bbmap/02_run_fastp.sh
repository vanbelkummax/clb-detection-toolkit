#!/bin/bash
set -e

# Usage: bash scripts/02_run_fastp.sh <SRR_ID>

SRR=$1
THREADS=2

INPUT_R1="data/raw_fastqs/${SRR}_R1.fastq.gz"
INPUT_R2="data/raw_fastqs/${SRR}_R2.fastq.gz"

OUTPUT_R1="pipeline1_read_mapping/fastp_qc/${SRR}_R1_trimmed.fastq.gz"
OUTPUT_R2="pipeline1_read_mapping/fastp_qc/${SRR}_R2_trimmed.fastq.gz"

HTML="pipeline1_read_mapping/fastp_qc/${SRR}_fastp.html"
JSON="pipeline1_read_mapping/fastp_qc/${SRR}_fastp.json"

LOG="logs/pipeline1/${SRR}_fastp.log"

echo "[$(date)] Running fastp for $SRR..." | tee "$LOG"

fastp \
  -i "$INPUT_R1" \
  -I "$INPUT_R2" \
  -o "$OUTPUT_R1" \
  -O "$OUTPUT_R2" \
  --detect_adapter_for_pe \
  --qualified_quality_phred 20 \
  --length_required 50 \
  --thread $THREADS \
  --html "$HTML" \
  --json "$JSON" \
  2>&1 | tee -a "$LOG"

echo "[$(date)] fastp complete for $SRR" | tee -a "$LOG"

# Create completion flag
touch "pipeline1_read_mapping/fastp_qc/${SRR}.fastp.done"
