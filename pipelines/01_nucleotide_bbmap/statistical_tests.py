#!/usr/bin/env python3
"""
Statistical tests for Pipeline 1 results:
1. Permutation test for CLB burden (mean Total RPM)
2. Fisher's exact test for CLB+ prevalence
"""

import pandas as pd
import numpy as np
from scipy import stats
from pathlib import Path

def permutation_test(group1, group2, n_permutations=100000):
    """Two-tailed permutation test for difference in means."""
    obs_mean1 = np.mean(group1)
    obs_mean2 = np.mean(group2)
    observed_diff = obs_mean1 - obs_mean2

    combined = np.concatenate([group1, group2])
    n1 = len(group1)

    perm_diffs = []
    for _ in range(n_permutations):
        np.random.shuffle(combined)
        perm_group1 = combined[:n1]
        perm_group2 = combined[n1:]
        perm_diff = np.mean(perm_group1) - np.mean(perm_group2)
        perm_diffs.append(perm_diff)

    perm_diffs = np.array(perm_diffs)
    p_value = np.mean(np.abs(perm_diffs) >= np.abs(observed_diff))

    # Effect size (Cohen's d)
    pooled_std = np.sqrt(((len(group1)-1)*np.var(group1, ddof=1) +
                          (len(group2)-1)*np.var(group2, ddof=1)) /
                         (len(group1) + len(group2) - 2))
    cohens_d = observed_diff / pooled_std if pooled_std > 0 else 0

    return {
        'observed_diff': observed_diff,
        'p_value': p_value,
        'cohens_d': cohens_d,
        'mean1': obs_mean1,
        'mean2': obs_mean2,
        'median1': np.median(group1),
        'median2': np.median(group2),
        'fold_change': obs_mean1 / obs_mean2 if obs_mean2 > 0 else np.inf
    }

def fishers_exact_test(neo_pos, neo_total, no_pos, no_total):
    """Fisher's exact test for 2×2 contingency table."""
    # Create contingency table
    # Rows: Neo-ABX, No-ABX
    # Cols: CLB+, CLB-
    table = [
        [neo_pos, neo_total - neo_pos],
        [no_pos, no_total - no_pos]
    ]

    oddsratio, p_value = stats.fisher_exact(table)

    return {
        'odds_ratio': oddsratio,
        'p_value': p_value,
        'neo_pos': neo_pos,
        'neo_total': neo_total,
        'neo_prevalence': neo_pos / neo_total * 100,
        'no_pos': no_pos,
        'no_total': no_total,
        'no_prevalence': no_pos / no_total * 100
    }

def main():
    # Load aggregated results
    summary_file = Path("pipeline1_read_mapping/summary_tables/pipeline1_sample_summary.tsv")
    df = pd.read_csv(summary_file, sep='\t')

    print("=" * 80)
    print("Pipeline 1 Statistical Analysis (bbmap, 95% identity)")
    print("=" * 80)
    print()

    # Split by exposure group
    neo_abx = df[df['Exposure'] == 'Neo-ABX']
    no_abx = df[df['Exposure'] == 'No-ABX']

    print(f"Sample sizes:")
    print(f"  Neo-ABX: n={len(neo_abx)}")
    print(f"  No-ABX:  n={len(no_abx)}")
    print()

    # ============================================================
    # 1. PERMUTATION TEST FOR CLB BURDEN (Total RPM)
    # ============================================================
    print("=" * 80)
    print("1. PERMUTATION TEST: CLB Burden (Total RPM)")
    print("=" * 80)
    print()

    neo_rpm = neo_abx['Total_CLB_RPM'].values
    no_rpm = no_abx['Total_CLB_RPM'].values

    perm_result = permutation_test(neo_rpm, no_rpm, n_permutations=100000)

    print(f"Neo-ABX:")
    print(f"  Mean:   {perm_result['mean1']:,.1f} RPM")
    print(f"  Median: {perm_result['median1']:,.1f} RPM")
    print()
    print(f"No-ABX:")
    print(f"  Mean:   {perm_result['mean2']:,.1f} RPM")
    print(f"  Median: {perm_result['median2']:,.1f} RPM")
    print()
    print(f"Difference: {perm_result['observed_diff']:+,.1f} RPM")
    print(f"Fold change: {perm_result['fold_change']:.2f}×")
    print(f"Cohen's d: {perm_result['cohens_d']:.3f}")
    print()
    print(f"Permutation test (100,000 iterations):")
    print(f"  P-value: {perm_result['p_value']:.6f}")

    if perm_result['p_value'] < 0.001:
        sig = "***"
    elif perm_result['p_value'] < 0.01:
        sig = "**"
    elif perm_result['p_value'] < 0.05:
        sig = "*"
    else:
        sig = "ns"
    print(f"  Significance: {sig}")
    print()

    # ============================================================
    # 2. FISHER'S EXACT TEST FOR CLB+ PREVALENCE
    # ============================================================
    print("=" * 80)
    print("2. FISHER'S EXACT TEST: CLB+ Prevalence (>9 genes with RPM>1)")
    print("=" * 80)
    print()

    neo_pos = (neo_abx['CLB_Status'] == 'Positive').sum()
    no_pos = (no_abx['CLB_Status'] == 'Positive').sum()

    fisher_result = fishers_exact_test(neo_pos, len(neo_abx), no_pos, len(no_abx))

    print(f"Neo-ABX: {fisher_result['neo_pos']}/{fisher_result['neo_total']} CLB+ ({fisher_result['neo_prevalence']:.1f}%)")
    print(f"No-ABX:  {fisher_result['no_pos']}/{fisher_result['no_total']} CLB+ ({fisher_result['no_prevalence']:.1f}%)")
    print()
    print(f"Fisher's exact test:")
    print(f"  Odds ratio: {fisher_result['odds_ratio']:.3f}")
    print(f"  P-value: {fisher_result['p_value']:.6f}")

    if fisher_result['p_value'] < 0.001:
        sig = "***"
    elif fisher_result['p_value'] < 0.01:
        sig = "**"
    elif fisher_result['p_value'] < 0.05:
        sig = "*"
    else:
        sig = "ns"
    print(f"  Significance: {sig}")
    print()

    # ============================================================
    # SAVE RESULTS
    # ============================================================
    results_df = pd.DataFrame([
        {
            'Test': 'Permutation (CLB Burden)',
            'Metric': 'Total RPM',
            'Neo_ABX_Mean': perm_result['mean1'],
            'No_ABX_Mean': perm_result['mean2'],
            'Difference': perm_result['observed_diff'],
            'Fold_Change': perm_result['fold_change'],
            'Effect_Size': perm_result['cohens_d'],
            'P_Value': perm_result['p_value'],
            'Significance': sig
        },
        {
            'Test': 'Fisher Exact (CLB+ Prevalence)',
            'Metric': 'CLB+ (>9 genes)',
            'Neo_ABX_Mean': fisher_result['neo_prevalence'],
            'No_ABX_Mean': fisher_result['no_prevalence'],
            'Difference': fisher_result['neo_prevalence'] - fisher_result['no_prevalence'],
            'Fold_Change': fisher_result['odds_ratio'],
            'Effect_Size': np.nan,
            'P_Value': fisher_result['p_value'],
            'Significance': sig
        }
    ])

    output_file = "pipeline1_read_mapping/summary_tables/pipeline1_statistical_tests.tsv"
    results_df.to_csv(output_file, sep='\t', index=False)
    print(f"Results saved to: {output_file}")
    print()

    # ============================================================
    # SUMMARY
    # ============================================================
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print()
    print("Pipeline 1 (bbmap, 95% nucleotide identity):")
    print(f"  • CLB burden {perm_result['fold_change']:.2f}× higher in Neo-ABX (p={perm_result['p_value']:.4f})")
    print(f"  • CLB+ prevalence {fisher_result['odds_ratio']:.2f}× higher odds in Neo-ABX (p={fisher_result['p_value']:.4f})")
    print()

if __name__ == '__main__':
    main()
