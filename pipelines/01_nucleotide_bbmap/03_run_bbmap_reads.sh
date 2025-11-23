#!/bin/bash
set -e

# Usage: bash scripts/03_run_bbmap_reads.sh <SRR_ID>

SRR=$1
THREADS=2

REF="data/reference/IHE3034_genome.fna"
INPUT_R1="pipeline1_read_mapping/fastp_qc/${SRR}_R1_trimmed.fastq.gz"
INPUT_R2="pipeline1_read_mapping/fastp_qc/${SRR}_R2_trimmed.fastq.gz"

OUTPUT_BAM="pipeline1_read_mapping/bbmap_alignments/${SRR}_aligned.bam"
COV_STATS="pipeline1_read_mapping/bbmap_alignments/${SRR}_coverage.txt"

LOG="logs/pipeline1/${SRR}_bbmap.log"

echo "[$(date)] Running bbmap for $SRR..." | tee "$LOG"

bbmap.sh \
  in="$INPUT_R1" \
  in2="$INPUT_R2" \
  ref="$REF" \
  out="$OUTPUT_BAM" \
  minid=0.95 \
  ambiguous=best \
  threads=$THREADS \
  covstats="$COV_STATS" \
  2>&1 | tee -a "$LOG"

# Sort and index BAM for samtools depth
echo "[$(date)] Sorting and indexing BAM..." | tee -a "$LOG"
samtools sort -@ $THREADS -o "pipeline1_read_mapping/bbmap_alignments/${SRR}_sorted.bam" "$OUTPUT_BAM"
samtools index "pipeline1_read_mapping/bbmap_alignments/${SRR}_sorted.bam"
mv "pipeline1_read_mapping/bbmap_alignments/${SRR}_sorted.bam" "$OUTPUT_BAM"
mv "pipeline1_read_mapping/bbmap_alignments/${SRR}_sorted.bam.bai" "pipeline1_read_mapping/bbmap_alignments/${SRR}_aligned.bam.bai"

echo "[$(date)] bbmap complete for $SRR" | tee -a "$LOG"

# Create completion flag
touch "pipeline1_read_mapping/bbmap_alignments/${SRR}.bbmap.done"
