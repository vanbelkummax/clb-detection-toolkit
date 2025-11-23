"""Generate publication-quality CLB island analysis figures at 95% identity.

Figures:
1. fig1_violin_final.(png|pdf) – log10 burden violin/box/strip plot by exposure.
2. fig2_prevalence_final.(png|pdf) – prevalence bars for key detection metrics.
3. fig6_dual_final.(png|pdf) – side-by-side prevalence for >9 genes vs complete island (standard).
4. fig5_burden_final.(png|pdf) – scatter of log10 burden vs island breadth, coloured by category.
5. fig_ecdf_final.(png|pdf) – ECDF of log10 burdens with median lines.
6. fig7_summary_final2.(png|pdf) – four-panel summary infographic.

The script expects the data files from clb_island_analysis_package/raw_data/ and
writes the figures to the current working directory. Run with:
    python generate_clb_figures.py
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")  # headless rendering

import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import fisher_exact
from statsmodels.distributions.empirical_distribution import ECDF
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SUMMARY_PATH = ROOT / "raw_data" / "clb_95_identity_summary.csv"
GENE_DETAILS_PATH = ROOT / "raw_data" / "clb_95_identity_gene_details.csv"


def load_data(summary_path: Path) -> pd.DataFrame:
    """Load sample-level summary and derive adjusted/log burdens."""
    summary = pd.read_csv(summary_path)
    pseudocount = 0.05  # avoids log(0)
    summary["adj_burden"] = summary["Normalized_Burden_per_10M"].replace(0, pseudocount)
    summary["log_burden"] = np.log10(summary["adj_burden"] + 1)
    return summary


def plot_violin(summary: pd.DataFrame, output_stem: str) -> None:
    counts = summary["Exposure_Group"].value_counts()
    mean_ratio = (
        summary.loc[summary["Exposure_Group"] == "Neo-ABX", "adj_burden"].mean()
        / summary.loc[summary["Exposure_Group"] == "No-ABX", "adj_burden"].mean()
    )
    ratio_text = f"{mean_ratio:.2f}× higher"
    palette = {"Neo-ABX": "#FF6B6B", "No-ABX": "#4ECDC4"}

    plt.figure(figsize=(6.5, 6), dpi=300)
    sns.violinplot(
        x="Exposure_Group",
        y="log_burden",
        data=summary,
        inner=None,
        palette=palette,
        cut=0,
        bw=0.2,
    )
    sns.boxplot(
        x="Exposure_Group",
        y="log_burden",
        data=summary,
        width=0.12,
        showcaps=True,
        boxprops={"facecolor": "white", "edgecolor": "black"},
        medianprops={"color": "black", "linewidth": 1.5},
        whiskerprops={"color": "black"},
        flierprops={"marker": "", "markersize": 0},
    )
    sns.stripplot(
        x="Exposure_Group",
        y="log_burden",
        data=summary,
        color="black",
        size=3,
        jitter=0.15,
        alpha=0.6,
    )

    max_y = summary["log_burden"].max()
    y_margin = 0.6
    y_bracket = max_y + y_margin * 0.3
    plt.plot(
        [0, 0, 1, 1],
        [y_bracket - 0.05, y_bracket, y_bracket, y_bracket - 0.05],
        color="black",
        linewidth=1.2,
    )
    plt.text(0.5, y_bracket + 0.05, "*", ha="center", va="bottom", fontsize=14)
    plt.text(0.5, y_bracket - 0.07, ratio_text, ha="center", va="top", fontsize=12)

    ax = plt.gca()
    xform = ax.get_xaxis_transform()
    ax.text(0, -0.12, f"n={counts['Neo-ABX']}", transform=xform, ha="center", va="top", fontsize=11, clip_on=False)
    ax.text(1, -0.12, f"n={counts['No-ABX']}", transform=xform, ha="center", va="top", fontsize=11, clip_on=False)

    plt.title("CLB Burden Distribution", fontsize=16, pad=15)
    plt.xlabel("")
    plt.ylabel("log$_{10}$(CLB burden + 1)", fontsize=14)
    plt.ylim(-0.3, max_y + y_margin)
    plt.yticks([tick for tick in plt.yticks()[0] if tick > 0], fontsize=12)
    plt.xticks([0, 1], ["Neo-ABX", "No-ABX"], fontsize=12)
    plt.tight_layout(rect=[0, 0.05, 1, 1])
    plt.savefig(f"{output_stem}.png", dpi=300)
    plt.savefig(f"{output_stem}.pdf")
    plt.close()


def _metric_counts(summary: pd.DataFrame, mask: pd.Series) -> tuple[int, int, int, int]:
    neo_mask = summary["Exposure_Group"] == "Neo-ABX"
    no_mask = summary["Exposure_Group"] == "No-ABX"
    neo_count = mask[neo_mask].sum()
    no_count = mask[no_mask].sum()
    return neo_count, int(neo_mask.sum()), no_count, int(no_mask.sum())


def plot_prevalence(summary: pd.DataFrame, output_stem: str) -> None:
    metrics = {
        ">9 genes": summary["Genes_With_Reads"] >= 10,
        "Complete island (≥99.5% breadth)": summary["Complete_Island_Standard"],
        "Burden ≥5k": summary["Complete_Island_Strict"],
        "Burden ≥5k": summary["Normalized_Burden_per_10M"] >= 5000,
        "Burden ≥50k": summary["Normalized_Burden_per_10M"] >= 50000,
    }

    fig, ax = plt.subplots(figsize=(7.5, 5), dpi=300)
    bar_width = 0.35
    x_positions = np.arange(len(metrics))

    for i, (label, mask) in enumerate(metrics.items()):
        neo_count, neo_total, no_count, no_total = _metric_counts(summary, mask)
        neo_prop = neo_count / neo_total
        no_prop = no_count / no_total

        ax.bar(i - bar_width / 2, neo_prop, width=bar_width, color="#FF6B6B")
        ax.bar(i + bar_width / 2, no_prop, width=bar_width, color="#4ECDC4")
        ax.text(i - bar_width / 2, neo_prop + 0.03, f"{neo_count}/{neo_total}", ha="center", fontsize=9)
        ax.text(i + bar_width / 2, no_prop + 0.03, f"{no_count}/{no_total}", ha="center", fontsize=9)

        contingency = np.array([[neo_count, neo_total - neo_count], [no_count, no_total - no_count]])
        _, p = fisher_exact(contingency)
        annotate = (p < 0.05) or (label == ">9 genes")
        if annotate:
            y_bar = max(neo_prop, no_prop) + 0.07
            ax.plot(
                [i - bar_width / 2, i - bar_width / 2, i + bar_width / 2, i + bar_width / 2],
                [y_bar - 0.01, y_bar, y_bar, y_bar - 0.01],
                color="black",
                linewidth=1,
            )
            stars = "**" if label == "Burden ≥50k" and p < 0.01 else "*"
            ax.text(i, y_bar + 0.02, stars, ha="center", va="bottom", fontsize=12)

    ax.set_xticks(x_positions)
    ax.set_xticklabels(list(metrics.keys()), rotation=15, ha="right", fontsize=11)
    ax.set_ylim(0, 0.5)
    ax.set_ylabel("Proportion of infants", fontsize=14)
    ax.set_title("CLB Detection Prevalence", fontsize=16)
    ax.legend(handles=[
        Patch(facecolor="#FF6B6B", edgecolor="#FF6B6B", label="Neo-ABX"),
        Patch(facecolor="#4ECDC4", edgecolor="#4ECDC4", label="No-ABX"),
    ], loc="upper right", fontsize=10)
    plt.tight_layout()
    plt.savefig(f"{output_stem}.png", dpi=300)
    plt.savefig(f"{output_stem}.pdf")
    plt.close()


def plot_dual(summary: pd.DataFrame, output_stem: str) -> None:
    metrics = {
        ">9 genes": summary["Genes_With_Reads"] >= 10,
        "Complete island (≥99.5% breadth)": summary["Complete_Island_Standard"],
    }

    fig, axes = plt.subplots(1, 2, figsize=(7, 4), dpi=300, sharey=True)
    for idx, (label, mask) in enumerate(metrics.items()):
        ax = axes[idx]
        neo_count, neo_total, no_count, no_total = _metric_counts(summary, mask)
        neo_prop = neo_count / neo_total
        no_prop = no_count / no_total
        ax.bar(0, neo_prop, width=0.4, color="#FF6B6B")
        ax.bar(1, no_prop, width=0.4, color="#4ECDC4")
        ax.text(0, neo_prop + 0.03, f"{neo_count}/{neo_total}", ha="center", fontsize=9)
        ax.text(1, no_prop + 0.03, f"{no_count}/{no_total}", ha="center", fontsize=9)
        contingency = np.array([[neo_count, neo_total - neo_count], [no_count, no_total - no_count]])
        _, p = fisher_exact(contingency)
        annotate = (p < 0.05) or (label == ">9 genes")
        if annotate:
            y_br = max(neo_prop, no_prop) + 0.08
            ax.plot([0, 0, 1, 1], [y_br - 0.01, y_br, y_br, y_br - 0.01], color="black", linewidth=1)
            ax.text(0.5, y_br + 0.02, "*", ha="center", va="bottom", fontsize=12)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Neo-ABX", "No-ABX"], fontsize=11)
        ax.set_title(label, fontsize=13)
    axes[0].set_ylabel("Proportion", fontsize=13)
    axes[0].set_ylim(0, 0.6)
    fig.suptitle("Comparison of Key Detection Criteria", fontsize=16)
    plt.tight_layout(rect=[0, 0.05, 1, 0.93])
    plt.savefig(f"{output_stem}.png", dpi=300)
    plt.savefig(f"{output_stem}.pdf")
    plt.close()


def plot_scatter(summary: pd.DataFrame, output_stem: str) -> None:
    summary = summary.copy()
    try:
        gene_df = pd.read_csv(GENE_DETAILS_PATH)
        breadth_sum = gene_df.groupby("SRR_ID")["Breadth_Coverage"].sum()
        summary["island_breadth"] = summary["SRR_ID"].map(breadth_sum).fillna(0)
    except FileNotFoundError:
        summary["island_breadth"] = 0.0

    def assign_category(row):
        if row["Normalized_Burden_per_10M"] >= 50000:
            return "Burden ≥50k"
        if row["Genes_With_Reads"] >= 10:
            return ">9 genes"
        if row["Normalized_Burden_per_10M"] >= 5000:
            return "Burden ≥5k"
        return "Negative"

    summary["category"] = summary.apply(assign_category, axis=1)

    palette = {
        "Burden ≥50k": "red",
        ">9 genes": "blue",
        "Burden ≥5k": "orange",
        "Negative": "gray",
    }

    fig, ax = plt.subplots(figsize=(7, 5), dpi=300)
    for category, group_df in summary.groupby("category"):
        ax.scatter(
            group_df["log_burden"],
            group_df["island_breadth"],
            label=category,
            color=palette[category],
            edgecolor="black",
            s=50,
        )
    ax.axhline(18, linestyle="--", color="black", alpha=0.6)
    xmin, _ = ax.get_xlim()
    ax.text(xmin, 18.0, "Complete island", va="bottom", ha="left", fontsize=9)
    ax.set_xlabel("log$_{10}$(CLB burden + 1)", fontsize=14)
    ax.set_ylabel("Total island breadth", fontsize=14)
    ax.set_title("CLB Burden vs. Island Breadth", fontsize=16)
    ax.legend(title="Category", fontsize=9, loc="lower right")
    plt.tight_layout()
    plt.savefig(f"{output_stem}.png", dpi=300)
    plt.savefig(f"{output_stem}.pdf")
    plt.close()


def plot_ecdf(summary: pd.DataFrame, output_stem: str) -> None:
    neo_vals = summary[summary["Exposure_Group"] == "Neo-ABX"]["log_burden"]
    no_vals = summary[summary["Exposure_Group"] == "No-ABX"]["log_burden"]
    ecdf_neo = ECDF(neo_vals)
    ecdf_no = ECDF(no_vals)
    median_neo = np.median(neo_vals)
    median_no = np.median(no_vals)

    fig, ax = plt.subplots(figsize=(6.5, 5), dpi=300)
    ax.step(ecdf_neo.x, ecdf_neo.y, where="post", label="Neo-ABX", color="#FF6B6B")
    ax.scatter(ecdf_neo.x, ecdf_neo.y, s=10, color="#FF6B6B", alpha=0.7)
    ax.step(ecdf_no.x, ecdf_no.y, where="post", label="No-ABX", color="#4ECDC4")
    ax.scatter(ecdf_no.x, ecdf_no.y, s=10, color="#4ECDC4", alpha=0.7)
    ax.axvline(median_neo, linestyle="--", color="#FF6B6B", label="Median Neo")
    ax.axvline(median_no, linestyle="--", color="#4ECDC4", label="Median No")
    ax.set_xlabel("log$_{10}$(CLB burden + 1)", fontsize=14)
    ax.set_ylabel("ECDF", fontsize=14)
    ax.set_title("Empirical Cumulative Distribution of CLB Burden", fontsize=16)
    ax.legend(fontsize=9, loc="lower right")
    plt.tight_layout()
    plt.savefig(f"{output_stem}.png", dpi=300)
    plt.savefig(f"{output_stem}.pdf")
    plt.close()


def plot_summary_infographic(summary: pd.DataFrame, output_stem: str) -> None:
    counts = summary["Exposure_Group"].value_counts()
    mean_ratio = (
        summary.loc[summary["Exposure_Group"] == "Neo-ABX", "adj_burden"].mean()
        / summary.loc[summary["Exposure_Group"] == "No-ABX", "adj_burden"].mean()
    )
    ratio_text = f"{mean_ratio:.2f}× higher"

    effect_data = [
        {"label": "Burden ≥50k", "OR": 7.42, "CI": (1.62, 34.00), "p": 0.0075},
        {"label": "Complete island", "OR": 3.43, "CI": (1.31, 8.98), "p": 0.0156},
        {"label": ">9 genes", "OR": 2.70, "CI": (1.14, 6.38), "p": 0.0355},
    ]

    fig = plt.figure(figsize=(10, 7), dpi=300)
    gs = gridspec.GridSpec(2, 2, figure=fig, height_ratios=[1, 1])
    palette = {"Neo-ABX": "#FF6B6B", "No-ABX": "#4ECDC4"}

    # Panel A
    ax1 = fig.add_subplot(gs[0, 0])
    sns.violinplot(ax=ax1, x="Exposure_Group", y="log_burden", data=summary, inner=None, palette=palette, cut=0, bw=0.2)
    sns.boxplot(
        ax=ax1,
        x="Exposure_Group",
        y="log_burden",
        data=summary,
        width=0.12,
        showcaps=True,
        boxprops={"facecolor": "white", "edgecolor": "black"},
        medianprops={"color": "black", "linewidth": 1.5},
        whiskerprops={"color": "black"},
        flierprops={"marker": "", "markersize": 0},
    )
    sns.stripplot(ax=ax1, x="Exposure_Group", y="log_burden", data=summary, color="black", size=2.5, jitter=0.1, alpha=0.5)
    max_y = summary["log_burden"].max()
    y_margin = 0.5
    y_bracket = max_y + y_margin * 0.25
    ax1.plot([0, 0, 1, 1], [y_bracket - 0.04, y_bracket, y_bracket, y_bracket - 0.04], color="black", linewidth=1)
    ax1.text(0.5, y_bracket + 0.04, "*", ha="center", va="bottom", fontsize=12)
    ax1.text(0.5, y_bracket - 0.05, ratio_text, ha="center", va="top", fontsize=10)
    ax1.set_title("Burden Distribution", fontsize=12)
    ax1.set_xlabel("")
    ax1.set_ylabel("log$_{10}$(CLB burden + 1)", fontsize=10)
    ax1.set_ylim(-0.35, max_y + y_margin)
    ax1.text(0, -0.3, f"n={counts['Neo-ABX']}", ha="center", va="center", fontsize=8)
    ax1.text(1, -0.3, f"n={counts['No-ABX']}", ha="center", va="center", fontsize=8)
    ax1.set_xticks([0, 1])
    ax1.set_xticklabels(["Neo-ABX", "No-ABX"], fontsize=9)
    ax1.tick_params(axis="y", labelsize=8)

    # Panel B
    ax2 = fig.add_subplot(gs[0, 1])
    key_metrics = {
        ">9 genes": summary["Genes_With_Reads"] >= 10,
        "Complete island": summary["Complete_Island_Standard"],
    }
    for i, (label, mask) in enumerate(key_metrics.items()):
        neo_count, neo_total, no_count, no_total = _metric_counts(summary, mask)
        neo_prop = neo_count / neo_total
        no_prop = no_count / no_total
        ax2.bar(i - 0.15, neo_prop, width=0.3, color="#FF6B6B")
        ax2.bar(i + 0.15, no_prop, width=0.3, color="#4ECDC4")
        ax2.text(i - 0.15, neo_prop + 0.02, f"{neo_count}/{neo_total}", ha="center", fontsize=6.5)
        ax2.text(i + 0.15, no_prop + 0.02, f"{no_count}/{no_total}", ha="center", fontsize=6.5)
        contingency = np.array([[neo_count, neo_total - neo_count], [no_count, no_total - no_count]])
        _, p = fisher_exact(contingency)
        if p < 0.05:
            y_br = max(neo_prop, no_prop) + 0.04
            ax2.plot([i - 0.15, i - 0.15, i + 0.15, i + 0.15], [y_br - 0.01, y_br, y_br, y_br - 0.01], color="black")
            ax2.text(i, y_br + 0.01, "*", ha="center", fontsize=10)
    ax2.set_xticks([0, 1])
    ax2.set_xticklabels(list(key_metrics.keys()), fontsize=8)
    ax2.set_ylim(0, 1.1)
    ax2.set_ylabel("Proportion", fontsize=8)
    ax2.set_title("Detection Prevalence", fontsize=10)

    # Panel C
    ax3 = fig.add_subplot(gs[1, 0])
    y_pos = np.arange(len(effect_data))
    ors = [d["OR"] for d in effect_data]
    ci_low = [d["CI"][0] for d in effect_data]
    ci_high = [d["CI"][1] for d in effect_data]
    ax3.errorbar(
        ors,
        y_pos,
        xerr=[np.array(ors) - np.array(ci_low), np.array(ci_high) - np.array(ors)],
        fmt="o",
        color="black",
        ecolor="black",
        capsize=3,
    )
    for i, d in enumerate(effect_data):
        star = "**" if d["p"] < 0.01 else "*" if d["p"] < 0.05 else ""
        ax3.text(d["OR"], i + 0.15, f"{d['OR']:.2f}{star}", ha="center", fontsize=8)
    ax3.axvline(1, linestyle="--", color="grey")
    ax3.set_yticks(y_pos)
    ax3.set_yticklabels([d["label"] for d in effect_data], fontsize=8)
    ax3.set_xscale("log")
    ax3.set_xlim(0.5, 40)
    ax3.set_xlabel("OR (log scale)", fontsize=8)
    ax3.set_title("Effect Sizes", fontsize=10)

    # Panel D
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.axis("off")
    neo_subset = summary[summary["Exposure_Group"] == "Neo-ABX"]
    no_subset = summary[summary["Exposure_Group"] == "No-ABX"]
    table = [
        ["Group", "n", "Median log burden", "Complete island", ">9 genes"],
        [
            "Neo-ABX",
            counts["Neo-ABX"],
            f"{np.median(neo_subset['log_burden']):.2f}",
            f"{int(neo_subset['Complete_Island_Standard'].sum())}/{counts['Neo-ABX']}",
            f"{int((neo_subset['Genes_With_Reads'] >= 10).sum())}/{counts['Neo-ABX']}",
        ],
        [
            "No-ABX",
            counts["No-ABX"],
            f"{np.median(no_subset['log_burden']):.2f}",
            f"{int(no_subset['Complete_Island_Standard'].sum())}/{counts['No-ABX']}",
            f"{int((no_subset['Genes_With_Reads'] >= 10).sum())}/{counts['No-ABX']}",
        ],
    ]
    for i, row in enumerate(table):
        y = 0.9 - i * 0.25
        for j, cell in enumerate(row):
            x = 0.05 + j * 0.22
            style = {"weight": "bold"} if i == 0 else {}
            ax4.text(x, y, str(cell), transform=ax4.transAxes, fontsize=7, **style)

    fig.suptitle("CLB Detection Summary", fontsize=14, y=0.95)
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    plt.savefig(f"{output_stem}.png", dpi=300)
    plt.savefig(f"{output_stem}.pdf")
    plt.close()


def main() -> None:
    summary = load_data(SUMMARY_PATH)
    plot_violin(summary, "fig1_violin_final")
    plot_prevalence(summary, "fig2_prevalence_final")
    plot_dual(summary, "fig6_dual_final")
    plot_scatter(summary, "fig5_burden_final")
    plot_ecdf(summary, "fig_ecdf_final")
    plot_summary_infographic(summary, "fig7_summary_final2")


if __name__ == "__main__":
    main()
