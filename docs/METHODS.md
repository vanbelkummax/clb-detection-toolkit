# Detailed Methodology

## Table of Contents
1. [Sample Collection & Metadata](#sample-collection--metadata)
2. [Pipeline 1: Nucleotide-Level Read Mapping](#pipeline-1-nucleotide-level-read-mapping)
3. [Pipeline 2: Protein-Level Detection](#pipeline-2-protein-level-detection)
4. [Pipeline 3: Assembly-Based MAG Reconstruction](#pipeline-3-assembly-based-mag-reconstruction)
5. [Statistical Analysis](#statistical-analysis)
6. [Computational Environment](#computational-environment)

---

## Sample Collection & Metadata

### Dataset Overview
- **Source Study:** Ryan et al. (2025) - Bifidobacteria support optimal infant vaccine responses
- **Sample Type:** Stool samples from 6-week-old infants
- **Sequencing:** Paired-end Illumina metagenomic shotgun sequencing
- **Total Samples:** 113 infants
  - Neo-ABX group: n=33 (neonatal antibiotic exposure)
  - No-ABX group: n=80 (no antibiotic exposure)

### Metadata Fields
```
- Sample ID (SRR accession)
- WCH ID (clinical identifier)
- Exposure group (Neo-ABX / No-ABX)
- Sequencing depth (total reads)
- Read length (typically 150bp paired-end)
```

---

## Pipeline 1: Nucleotide-Level Read Mapping

### Overview
**Purpose:** Fast prevalence screening using direct nucleotide alignment
**Algorithm:** BBMap v39.37
**Runtime:** ~2 hours for 113 samples (parallel processing)
**Best for:** Initial screening, prevalence estimation

### Step-by-Step Protocol

#### 1. Quality Control (fastp v0.23.4)
```bash
fastp \
    --in1 ${SAMPLE}_1.fastq.gz \
    --in2 ${SAMPLE}_2.fastq.gz \
    --out1 ${SAMPLE}_clean_1.fastq.gz \
    --out2 ${SAMPLE}_clean_2.fastq.gz \
    --qualified_quality_phred 30 \
    --unqualified_percent_limit 40 \
    --n_base_limit 5 \
    --length_required 50 \
    --thread 4
```

**Parameters:**
- Quality threshold: Q30
- Minimum length: 50bp
- Maximum N bases: 5
- Adapter auto-detection enabled

#### 2. CLB Quantification (bbmap v39.37)
```bash
bbmap.sh \
    in1=${SAMPLE}_clean_1.fastq.gz \
    in2=${SAMPLE}_clean_2.fastq.gz \
    ref=IHE3034_clb_island.fasta \
    outm=${SAMPLE}_clb_mapped.sam \
    covstats=${SAMPLE}_clb_coverage.txt \
    rpkm=${SAMPLE}_clb_rpkm.txt \
    minid=0.95 \
    ambiguous=toss \
    threads=10
```

**Parameters:**
- `minid=0.95`: 95% nucleotide identity threshold
- `ambiguous=toss`: Discard ambiguously mapped reads
- Reference: IHE3034 complete *pks* island (5,108,383 bp)

#### 3. RPM Calculation
```python
# Extract reads per gene
for gene in CLB_GENES:
    gene_reads = count_mapped_reads(gene)
    total_microbial_reads = get_total_reads() - human_reads
    rpm = (gene_reads / total_microbial_reads) * 1_000_000
```

**Normalization:**
- RPM = Reads Per Million microbial reads
- Accounts for sequencing depth variation
- Human reads excluded from denominator

#### 4. Gene Detection Criteria
- **CLB+:** >9 genes with RPM > 1
- Rationale: Conservative threshold minimizing false positives

---

## Pipeline 2: Protein-Level Detection

### Overview
**Purpose:** High-sensitivity functional annotation with multi-threshold validation
**Algorithm:** DIAMOND blastx v2.1.10
**Innovation:** One-pass search + post-hoc filtering (4× faster than multiple runs)
**Runtime:** ~6 hours for 113 samples

### Step-by-Step Protocol

#### 1. Quality Trimming (BBDuk v39.10)
```bash
bbduk.sh \
    in1=${SAMPLE}_1.fastq.gz \
    in2=${SAMPLE}_2.fastq.gz \
    out1=${SAMPLE}_trimmed_1.fastq.gz \
    out2=${SAMPLE}_trimmed_2.fastq.gz \
    ref=adapters.fa \
    ktrim=r k=23 mink=11 hdist=1 \
    qtrim=rl trimq=20 \
    minlen=30 \
    tpe tbo
```

**Parameters:**
- Adapter trimming: k-mer based (k=23)
- Quality trimming: Q20 threshold
- Minimum length: 30bp
- Paired-end awareness enabled

#### 2. Human DNA Depletion (Bowtie2 v2.5.4)
```bash
bowtie2 \
    --very-sensitive-local \
    -x GRCh38 \
    -1 ${SAMPLE}_trimmed_1.fastq.gz \
    -2 ${SAMPLE}_trimmed_2.fastq.gz \
    -S ${SAMPLE}_human.sam \
    --un-conc-gz ${SAMPLE}_microbial_%.fastq.gz \
    --threads 10
```

**Parameters:**
- Reference: GRCh38 (human genome)
- Mode: Very-sensitive-local
- Output: Unmapped reads (microbial fraction)
- Result: <1% human contamination on average

#### 3. CLB Gene Detection (DIAMOND blastx)
```bash
diamond blastx \
    --query ${SAMPLE}_microbial_1.fastq.gz \
    --db clb_island_complete.dmnd \
    --out ${SAMPLE}_clb_hits.m8 \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
    --id 70 \
    --evalue 1e-5 \
    --max-target-seqs 25 \
    --fast \
    --block-size 6.0 \
    --threads 10
```

**Parameters:**
- Database: 18 CLB proteins (clbA-clbR)
- Initial threshold: 70% amino acid identity
- E-value cutoff: 1e-5
- Max targets per query: 25
- Fast mode enabled for speed

**Database Construction:**
```bash
diamond makedb \
    --in clb_island_complete.fasta \
    --db clb_island_complete
```

#### 4. Multi-Threshold Post-Hoc Filtering
```python
# Single DIAMOND run, filter to multiple thresholds
thresholds = [70, 75, 80, 85, 90, 95, 100]

for threshold in thresholds:
    filtered_hits = hits[hits['pident'] >= threshold]
    # Calculate metrics at each threshold
    unique_reads = count_unique_queries(filtered_hits)
    total_reads = count_all_alignments(filtered_hits)
    breadth = calculate_nucleotide_breadth(filtered_hits)
```

**Efficiency Gain:** 4× faster than running DIAMOND 7 times

#### 5. Breadth Calculation (Nucleotide-Level)
```python
def calculate_breadth(alignments, gene_length):
    """
    Calculate breadth as fraction of gene covered by alignments
    at nucleotide level (not amino acid level)
    """
    covered_positions = set()
    for aln in alignments:
        # Convert protein positions to nucleotide
        nt_start = aln.sstart * 3
        nt_end = aln.send * 3
        covered_positions.update(range(nt_start, nt_end))

    breadth_fraction = len(covered_positions) / (gene_length * 3)
    return breadth_fraction
```

#### 6. Gene Detection Criteria
- **Minimum unique reads:** ≥5 per gene
- **Minimum breadth:** ≥10% (nucleotide fraction)
- **Short genes** (clbA, clbD, clbE, clbQ, clbR): ≥20aa minimum alignment
- **Other genes:** ≥30aa minimum alignment

**Complete Island Definitions:**
1. **Standard:** All 18 genes with ≥5 reads & ≥10% breadth
2. **Strict:** All 18 genes with ≥99.9% breadth

#### 7. Normalization
```python
burden_per_10M = (total_clb_reads / total_microbial_reads) * 10_000_000
```

---

## Pipeline 3: Assembly-Based MAG Reconstruction

### Overview
**Purpose:** Prove genome-coherent co-location of *pks* + virulence genes
**Approach:** Assembly → Binning → QC → Annotation
**Runtime:** 10-14 days for 398 assemblies

### Step-by-Step Protocol

#### 1. Metagenomic Assembly (metaSPAdes v3.15+)
```bash
metaspades.py \
    -1 ${SAMPLE}_1.fastq.gz \
    -2 ${SAMPLE}_2.fastq.gz \
    -o ${SAMPLE}_assembly \
    --threads 16 \
    --memory 250
```

**Parameters:**
- K-mer sizes: Auto-detected (typically 21, 33, 55, 77)
- Memory: 250 GB
- Output: Scaffolds ≥500bp

**Input:** 398 metaSPAdes assemblies from Feargal Ryan

#### 2. MAG Binning (Multi-Tool Consensus)

**MetaBAT2 v2.15:**
```bash
metabat2 \
    -i ${SAMPLE}.scaffolds.fasta \
    -o bins_metabat2/${SAMPLE} \
    --minContig 2500 \
    --minCV 1.0 \
    --minCVSum 1.0
```

**MaxBin2 v2.2.7:**
```bash
run_MaxBin.pl \
    -contig ${SAMPLE}.scaffolds.fasta \
    -out bins_maxbin2/${SAMPLE} \
    -min_contig_length 2500 \
    -thread 16
```

**SemiBin2 v2.0:**
```bash
SemiBin2 single_easy_bin \
    -i ${SAMPLE}.scaffolds.fasta \
    -o bins_semibin2 \
    --environment human_gut \
    -p 16
```

**DAS_Tool v1.1.5 (Consensus):**
```bash
DAS_Tool \
    -i metabat2_scaffolds2bin.tsv,maxbin2_scaffolds2bin.tsv,semibin2_scaffolds2bin.tsv \
    -l MetaBAT2,MaxBin2,SemiBin2 \
    -c ${SAMPLE}.scaffolds.fasta \
    -o ${SAMPLE}_dastool \
    --write_bins \
    --threads 16
```

#### 3. Quality Control

**CheckM2 v1.0+ (Completeness/Contamination):**
```bash
checkm2 predict \
    --input bins/ \
    --output-directory checkm2_results \
    --threads 16 \
    --extension .fa
```

**Quality Thresholds:**
- High-quality: Completeness ≥90%, Contamination ≤5%
- Medium-quality: Completeness ≥50%, Contamination ≤10%

**GUNC v1.0.5 (Chimera Detection):**
```bash
gunc run \
    --input_dir bins/ \
    --out_dir gunc_results \
    --threads 16 \
    --db_file ~/.gunc/gunc_db_progenomes2.1.dmnd
```

**Chimera Threshold:** GUNC score > 0.45 flagged for review

#### 4. Taxonomic Classification (GTDB-Tk R220)
```bash
gtdbtk classify_wf \
    --genome_dir bins/ \
    --out_dir gtdb_results \
    --extension .fa \
    --cpus 16 \
    --pplacer_cpus 4
```

**Database:** GTDB R220 (release 220)

#### 5. Gene Calling (Prodigal v2.6.3)
```bash
prodigal \
    -i ${MAG}.fa \
    -a ${MAG}.faa \
    -d ${MAG}.fna \
    -f gff \
    -p meta
```

**Mode:** Metagenomic (`-p meta`)

#### 6. CLB Island Detection (DIAMOND blastp)
```bash
diamond blastp \
    --query ${MAG}.faa \
    --db clb_island_proteins.dmnd \
    --out ${MAG}_clb_hits.tsv \
    --outfmt 6 \
    --id 95 \
    --evalue 1e-10 \
    --max-target-seqs 1 \
    --threads 4
```

#### 7. Synteny Verification
```python
def check_synteny(mag_clb_genes):
    """
    Verify clbB-P-Q order (intact pks island marker)
    """
    required_order = ['clbB', 'clbC', 'clbD', ..., 'clbP', 'clbQ']

    # Check if genes appear in correct order on scaffolds
    for scaffold in mag.scaffolds:
        genes_on_scaffold = sorted(clb_genes, key=lambda x: x.position)
        if genes_match_order(genes_on_scaffold, required_order):
            return True
    return False
```

#### 8. Virulence Toolkit Annotation
**64 additional virulence/fitness genes:**
- Type 1 fimbriae (fim operon): 7 genes
- E. coli common pilus (ecp): 6 genes
- Respiration (napA, narG, cyd, cyo): 8 genes
- β-oxidation (fad genes): 6 genes
- K1 capsule: 15 genes
- Adhesins (sfa/foc, dr/afa, ompA): 8 genes

---

## Statistical Analysis

### Primary Tests

#### 1. Permutation Testing (CLB Burden)
```python
import numpy as np
from scipy import stats

def permutation_test(group1, group2, n_permutations=100000, seed=42):
    """
    Non-parametric permutation test for difference in means
    """
    np.random.seed(seed)

    # Observed difference
    observed_diff = np.mean(group1) - np.mean(group2)

    # Combine groups
    combined = np.concatenate([group1, group2])
    n1 = len(group1)

    # Permutation distribution
    perm_diffs = []
    for _ in range(n_permutations):
        np.random.shuffle(combined)
        perm_group1 = combined[:n1]
        perm_group2 = combined[n1:]
        perm_diff = np.mean(perm_group1) - np.mean(perm_group2)
        perm_diffs.append(perm_diff)

    # Two-tailed p-value
    perm_diffs = np.array(perm_diffs)
    p_value = np.mean(np.abs(perm_diffs) >= np.abs(observed_diff))

    return p_value, observed_diff
```

**Parameters:**
- Iterations: 100,000
- Random seed: 42 (reproducibility)
- Test: Two-tailed

#### 2. Fisher's Exact Test (Prevalence)
```R
# With Haldane-Anscombe correction
fisher_test_corrected <- function(a, b, c, d) {
    # Add 0.5 to all cells
    contingency_table <- matrix(c(a+0.5, b+0.5, c+0.5, d+0.5), nrow=2)
    test <- fisher.test(contingency_table)
    return(test)
}
```

**Correction:** Haldane-Anscombe (add 0.5 to all cells)
**Rationale:** Improves performance with small counts

#### 3. Effect Sizes

**Cohen's d:**
```python
def cohens_d(group1, group2):
    """
    Standardized effect size
    """
    mean1, mean2 = np.mean(group1), np.mean(group2)
    std1, std2 = np.std(group1, ddof=1), np.std(group2, ddof=1)
    n1, n2 = len(group1), len(group2)

    # Pooled standard deviation
    pooled_std = np.sqrt(((n1-1)*std1**2 + (n2-1)*std2**2) / (n1 + n2 - 2))

    d = (mean1 - mean2) / pooled_std
    return d
```

**Odds Ratio with 95% CI:**
```python
from statsmodels.stats.contingency_tables import Table2x2

def odds_ratio_ci(a, b, c, d):
    """
    OR with 95% confidence interval
    """
    table = Table2x2([[a, b], [c, d]])
    or_value = table.oddsratio
    ci_low, ci_high = table.oddsratio_confint()
    return or_value, (ci_low, ci_high)
```

#### 4. Multiple Testing Correction
```python
from statsmodels.stats.multitest import multipletests

def fdr_correction(pvalues, alpha=0.05):
    """
    Benjamini-Hochberg FDR correction
    """
    reject, pvals_corrected, _, _ = multipletests(
        pvalues,
        alpha=alpha,
        method='fdr_bh'
    )
    return pvals_corrected, reject
```

**Method:** Benjamini-Hochberg FDR
**Threshold:** α = 0.05

---

## Computational Environment

### Hardware Requirements

**Minimum Specifications:**
- **RAM:** 512 GB (for MAG binning)
- **Storage:** 8 TB
  - Raw data: ~4 TB
  - Databases: 81.5 GB (GUNC + GTDB-Tk)
  - Results: ~3 TB
- **CPU:** 16+ cores recommended
- **OS:** Linux (Ubuntu 20.04+ or CentOS 7+)

**Recommended:**
- SLURM cluster with checkpointed execution
- High-performance computing (HPC) environment

### Software Versions

**Alignment & Quantification:**
- bbmap v39.37
- DIAMOND v2.1.10
- Bowtie2 v2.5.4
- fastp v0.23.4

**Assembly & Binning:**
- metaSPAdes v3.15+
- MetaBAT2 v2.15
- MaxBin2 v2.2.7
- SemiBin2 v2.0
- DAS_Tool v1.1.5

**Quality Control:**
- CheckM2 v1.0+
- GUNC v1.0.5
- GTDB-Tk R220

**Gene Calling & Annotation:**
- Prodigal v2.6.3

**Statistical Analysis:**
- Python 3.8+ (numpy 1.24+, scipy 1.10+, pandas 2.0+, statsmodels 0.14+)
- R 4.3+ (tidyverse 2.0+, logistf 1.25+)

**Visualization:**
- matplotlib v3.7+
- seaborn v0.12+

### Conda Environments

See `envs/` directory for complete environment specifications:
- `pipeline1_bbmap.yml`
- `pipeline2_diamond.yml`
- `pipeline3_assembly.yml`
- `analysis_stats.yml`

---

## Quality Control Metrics

### Sequencing Quality
- **Read count:** 5M - 50M reads per sample
- **Read length:** 150bp paired-end
- **Q30 score:** >80% of bases
- **Human contamination:** <1% (after depletion)

### MAG Quality
- **High-quality MAGs:** Completeness ≥90%, Contamination ≤5%
- **Medium-quality MAGs:** Completeness ≥50%, Contamination ≤10%
- **Chimera detection:** GUNC score ≤0.45

### Pipeline Validation
- **Positive controls:** E. coli IHE3034 (*pks+* reference)
- **Negative controls:** Human DNA-only samples
- **Technical replicates:** ≥3 samples processed in duplicate

---

## Reproducibility

**All analyses are fully reproducible:**
1. Conda environments lock software versions
2. Random seeds specified for stochastic processes
3. Reference databases versioned (GTDB R220, GRCh38)
4. All parameters documented
5. Scripts available in this repository

---

## References

1. Kang DD, et al. (2019) MetaBAT 2: an adaptive binning algorithm for robust and efficient genome reconstruction from metagenome assemblies. *PeerJ* 7:e7359.

2. Wu YW, et al. (2016) MaxBin 2.0: an automated binning algorithm to recover genomes from multiple metagenomic datasets. *Bioinformatics* 32:605-607.

3. Pan S, et al. (2022) SemiBin2: a semi-supervised contigs binning method for recovering high-quality metagenome-assembled genomes. *Genome Biology* 23:102.

4. Sieber CMK, et al. (2018) Recovery of genomes from metagenomes via a dereplication, aggregation and scoring strategy. *Nature Microbiology* 3:836-843.

5. Chklovski A, et al. (2023) CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning. *Nature Methods* 20:1203-1212.

6. Orakov A, et al. (2021) GUNC: detection of chimerism and contamination in prokaryotic genomes. *Genome Biology* 22:178.

7. Chaumeil PA, et al. (2022) GTDB-Tk v2: memory friendly classification with the Genome Taxonomy Database. *Bioinformatics* 38:5315-5316.

8. Buchfink B, et al. (2021) Sensitive protein alignments at tree-of-life scale using DIAMOND. *Nature Methods* 18:366-368.

9. Langmead B & Salzberg SL (2012) Fast gapped-read alignment with Bowtie 2. *Nature Methods* 9:357-359.

10. Hyatt D, et al. (2010) Prodigal: prokaryotic gene recognition and translation initiation site identification. *BMC Bioinformatics* 11:119.

---

*For questions about methodology, please open an issue on GitHub.*
