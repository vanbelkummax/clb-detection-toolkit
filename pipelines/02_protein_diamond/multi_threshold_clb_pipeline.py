#!/usr/bin/env python3
"""
Multi-threshold CLB Island Analysis Pipeline v2.0
Fixed version with correct database, breadth calculation, and dual complete island metrics
"""

import os
import sys
import json
import subprocess
import gzip
from datetime import datetime
from collections import defaultdict
import argparse

# Pipeline configuration
BBDUK_PATH = "/mnt/x/validation_clbB/bbmap/bbduk.sh"
HUMAN_DB = "/mnt/x/validation_clbB/databases/GRCh38"
CLB_DB = "/mnt/x/validation_clbB/proteins/clb_island_complete.dmnd"  # FIXED: Use full 18-gene database
DIAMOND_PATH = "/mnt/x/validation_clbB/bin/diamond"
ADAPTERS = "/mnt/x/validation_clbB/bbmap/resources/adapters.fa"

# CLB gene metadata (nucleotide lengths)
CLB_GENES = {
    'clbA': {'length': 1521, 'short': True},
    'clbB': {'length': 23817, 'short': False},
    'clbC': {'length': 7503, 'short': False},
    'clbD': {'length': 1242, 'short': True},
    'clbE': {'length': 726, 'short': True},
    'clbF': {'length': 2319, 'short': False},
    'clbG': {'length': 2814, 'short': False},
    'clbH': {'length': 9408, 'short': False},
    'clbI': {'length': 6060, 'short': False},
    'clbJ': {'length': 12522, 'short': False},
    'clbK': {'length': 13017, 'short': False},
    'clbL': {'length': 3846, 'short': False},
    'clbM': {'length': 3150, 'short': False},
    'clbN': {'length': 9798, 'short': False},
    'clbO': {'length': 5343, 'short': False},
    'clbP': {'length': 3012, 'short': False},
    'clbQ': {'length': 1215, 'short': True},
    'clbR': {'length': 594, 'short': True}
}

def run_command(cmd, description="Running command"):
    """Execute shell command and handle errors."""
    print(f"{description}...")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error in {description}:")
        print(result.stderr)
        return None
    return result.stdout

def quality_trim(r1_path, r2_path, output_dir, threads=10):
    """Run BBDuk for quality trimming and adapter removal."""
    trimmed_r1 = os.path.join(output_dir, "trimmed_R1.fastq.gz")
    trimmed_r2 = os.path.join(output_dir, "trimmed_R2.fastq.gz")
    stats_file = os.path.join(output_dir, "bbduk_stats.txt")

    cmd = (
        f"{BBDUK_PATH} in={r1_path} in2={r2_path} "
        f"out={trimmed_r1} out2={trimmed_r2} "
        f"ref={ADAPTERS} "
        f"ktrim=r k=23 mink=11 hdist=1 "
        f"qtrim=rl trimq=20 minlen=30 "
        f"tpe tbo threads={threads} "
        f"stats={stats_file} 2>/dev/null"
    )

    if run_command(cmd, "Quality trimming with BBDuk") is None:
        return None, None

    return trimmed_r1, trimmed_r2

def human_depletion(r1_path, r2_path, output_dir, sample_name, threads=10):
    """Remove human reads using Bowtie2."""
    # Check for appropriate index files
    if not (os.path.exists(f"{HUMAN_DB}.1.bt2") or os.path.exists(f"{HUMAN_DB}.1.bt2l")):
        print(f"Error: Bowtie2 index not found at {HUMAN_DB}")
        return None, None

    depleted_r1 = os.path.join(output_dir, "microbial_R1.fastq.gz")
    depleted_r2 = os.path.join(output_dir, "microbial_R2.fastq.gz")

    cmd = (
        f"bowtie2 -x {HUMAN_DB} -1 {r1_path} -2 {r2_path} "
        f"--threads {threads} --very-sensitive-local -S /dev/null 2>/dev/null "
        f"| samtools view -f 4 -F 256 - "
        f"| samtools fastq -1 {depleted_r1} -2 {depleted_r2} -s /dev/null -"
    )

    if run_command(cmd, "Human read depletion with Bowtie2") is None:
        return None, None

    # Count microbial reads
    count_cmd = f"zcat {depleted_r1} | wc -l"
    result = subprocess.run(count_cmd, shell=True, capture_output=True, text=True)
    microbial_reads = int(result.stdout.strip()) // 4

    return depleted_r1, depleted_r2, microbial_reads

def run_diamond(r1_path, r2_path, output_dir, threads=10):
    """Run DIAMOND blastx against CLB database."""
    # Combine paired-end reads for DIAMOND
    combined_fastq = os.path.join(output_dir, "combined.fastq.gz")
    cmd = f"cat {r1_path} {r2_path} > {combined_fastq}"
    if run_command(cmd, "Combining FASTQ files") is None:
        return None

    # Run DIAMOND at 70% identity (will post-slice to other thresholds)
    diamond_output = os.path.join(output_dir, "clb_hits_id70_raw.m8")
    diamond_cmd = (
        f"{DIAMOND_PATH} blastx --db {CLB_DB} --query {combined_fastq} "
        f"--out {diamond_output} --threads {threads} "
        f"--id 70 --evalue 1e-5 --max-target-seqs 25 "
        f"--block-size 6.0 --fast --outfmt 6 qseqid sseqid pident length "
        f"mismatch gapopen qstart qend sstart send evalue bitscore qlen"
    )

    if run_command(diamond_cmd, "DIAMOND search") is None:
        return None

    # Clean up combined file
    os.remove(combined_fastq)

    return diamond_output

def slice_bins_and_count(m8_path, sample_dir, thresholds=(70, 75, 80, 85, 90, 95, 100),
                        min_alignment_short=20, min_alignment_long=30,
                        min_reads=5, min_breadth=0.10):  # FIXED: breadth as fraction
    """Slice DIAMOND results by identity thresholds and count gene detection."""

    threshold_stats = {}

    for threshold in thresholds:
        # Read and filter hits by identity threshold
        gene_hits = defaultdict(lambda: {
            'unique_reads': set(),
            'total_reads': defaultdict(int),  # FIXED: Track total separately
            'positions': set()
        })

        with open(m8_path) as f:
            for line in f:
                fields = line.strip().split('\t')
                read_id = fields[0]
                gene_id = fields[1].split('_')[0]  # Extract gene name
                pident = float(fields[2])
                align_len = int(fields[3])
                sstart = int(fields[8])
                send = int(fields[9])

                # Skip if below threshold
                if pident < threshold:
                    continue

                # Check minimum alignment length
                if gene_id in CLB_GENES:
                    min_len = min_alignment_short if CLB_GENES[gene_id]['short'] else min_alignment_long
                    if align_len < min_len:
                        continue
                else:
                    continue  # Skip unknown genes

                # Record hit
                gene_hits[gene_id]['unique_reads'].add(read_id)
                gene_hits[gene_id]['total_reads'][read_id] += 1  # Count multiple alignments

                # Record positions for breadth calculation
                # Convert protein positions to nucleotide positions (multiply by 3)
                nt_start = min(sstart, send) * 3 - 2
                nt_end = max(sstart, send) * 3
                for pos in range(nt_start, nt_end + 1):
                    gene_hits[gene_id]['positions'].add(pos)

        # Calculate statistics for this threshold
        stats = {
            'identity_threshold': threshold,
            'genes': {},
            'genes_detected': 0,
            'complete_island': False,
            'complete_island_strict': False,  # NEW: ≥99.9% breadth for all genes
            'unique_reads': 0,
            'total_pairs': 0
        }

        # Process each gene
        genes_with_high_breadth = 0
        for gene_name, gene_info in CLB_GENES.items():
            gene_data = gene_hits.get(gene_name, {
                'unique_reads': set(),
                'total_reads': defaultdict(int),
                'positions': set()
            })

            unique = len(gene_data['unique_reads'])
            total = sum(gene_data['total_reads'].values())

            # Calculate breadth as fraction (0-1)
            breadth = len(gene_data['positions']) / gene_info['length'] if gene_info['length'] > 0 else 0

            # Detection based on reads AND breadth
            detected = unique >= min_reads and breadth >= min_breadth

            stats['genes'][gene_name] = {
                'unique_reads': unique,
                'total_reads': total,
                'breadth': breadth,  # FIXED: Store as fraction not percentage
                'mean_depth': total / gene_info['length'] if gene_info['length'] > 0 else 0,
                'detected': detected
            }

            if detected:
                stats['genes_detected'] += 1

            # Check for high breadth (≥99.9%)
            if breadth >= 0.999:
                genes_with_high_breadth += 1

            stats['unique_reads'] += unique
            stats['total_pairs'] += total

        # Check if complete island (standard: all 18 genes detected)
        stats['complete_island'] = stats['genes_detected'] == 18

        # Check if complete island (strict: all 18 genes with ≥99.9% breadth)
        stats['complete_island_strict'] = genes_with_high_breadth == 18

        threshold_stats[f'id{threshold}'] = stats

    return threshold_stats

def process_sample(sample_name, r1_path, r2_path, output_dir, exposure_group=None):
    """Process a single sample through the complete pipeline."""

    print(f"\n{'='*60}")
    print(f"Processing sample: {sample_name}")
    print(f"{'='*60}")

    sample_dir = os.path.join(output_dir, sample_name)
    os.makedirs(sample_dir, exist_ok=True)

    # Step 1: Quality trimming
    trimmed_r1, trimmed_r2 = quality_trim(r1_path, r2_path, sample_dir)
    if not trimmed_r1:
        return None

    # Step 2: Human depletion
    depleted_r1, depleted_r2, microbial_reads = human_depletion(
        trimmed_r1, trimmed_r2, sample_dir, sample_name
    )
    if not depleted_r1:
        return None

    print(f"Microbial reads: {microbial_reads:,}")

    # Step 3: DIAMOND search
    diamond_output = run_diamond(depleted_r1, depleted_r2, sample_dir)
    if not diamond_output:
        return None

    # Step 4: Multi-threshold analysis
    threshold_stats = slice_bins_and_count(
        diamond_output, sample_dir,
        thresholds=(70, 75, 80, 85, 90, 95, 100)
    )

    # Calculate normalized burden for each threshold
    for threshold_data in threshold_stats.values():
        if microbial_reads > 0:
            threshold_data['normalized_burden'] = (
                threshold_data['unique_reads'] / microbial_reads * 10000000
            )
        else:
            threshold_data['normalized_burden'] = 0

    # Create summary JSON
    summary = {
        'sample': sample_name,
        'exposure': exposure_group,
        'microbial_reads': microbial_reads,
        'processing_date': datetime.now().isoformat(),
        'pipeline_version': 'multi-threshold-2.0',
        'thresholds': threshold_stats
    }

    # Save summary
    stats_dir = os.path.join(output_dir, 'stats')
    os.makedirs(stats_dir, exist_ok=True)

    summary_path = os.path.join(stats_dir, f'{sample_name}_summary.json')
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)

    print(f"Summary saved to: {summary_path}")

    # Print 95% identity results
    if 'id95' in threshold_stats:
        id95 = threshold_stats['id95']
        print(f"\n95% Identity Results:")
        print(f"  Genes detected: {id95['genes_detected']}/18")
        print(f"  Complete island (standard): {id95['complete_island']}")
        print(f"  Complete island (≥99.9% breadth): {id95['complete_island_strict']}")
        print(f"  Normalized burden: {id95['normalized_burden']:.1f} reads/10M")

    # Clean up intermediate files to save space
    for f in [trimmed_r1, trimmed_r2, depleted_r1, depleted_r2]:
        if os.path.exists(f):
            os.remove(f)

    return summary

def main():
    parser = argparse.ArgumentParser(
        description='Multi-threshold CLB Island Analysis Pipeline v2.0',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Pipeline performs:
1. Quality trimming with BBDuk (Q20, minlen=30, adapter removal)
2. Human read depletion with Bowtie2 (--very-sensitive-local)
3. DIAMOND BLASTx search against full CLB island (18 genes, 70% identity)
4. Post-slicing to multiple thresholds: 70, 75, 80, 85, 90, 95, 100%
5. Gene detection and burden calculation
6. Dual complete island metrics (standard and ≥99.9% breadth)

Example:
    python multi_threshold_clb_pipeline_v2.py SRR12345678 \\
        /data/SRR12345678_1.fastq.gz /data/SRR12345678_2.fastq.gz \\
        /output/results --exposure Neo-ABX
        """
    )

    parser.add_argument('sample_name', help='Sample identifier (e.g., SRR12345678)')
    parser.add_argument('r1_fastq', help='Path to R1 FASTQ file (can be gzipped)')
    parser.add_argument('r2_fastq', help='Path to R2 FASTQ file (can be gzipped)')
    parser.add_argument('output_dir', help='Output directory for results')
    parser.add_argument('--exposure', choices=['Neo-ABX', 'No-ABX'],
                       help='Exposure group for the sample')

    args = parser.parse_args()

    # Validate input files
    for f in [args.r1_fastq, args.r2_fastq]:
        if not os.path.exists(f):
            print(f"Error: Input file not found: {f}")
            sys.exit(1)

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Process sample
    result = process_sample(
        args.sample_name,
        args.r1_fastq,
        args.r2_fastq,
        args.output_dir,
        args.exposure
    )

    if result:
        print("\n✓ Pipeline completed successfully!")
    else:
        print("\n✗ Pipeline failed!")
        sys.exit(1)

if __name__ == '__main__':
    main()