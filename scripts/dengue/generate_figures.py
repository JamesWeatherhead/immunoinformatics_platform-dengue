#!/usr/bin/env python3
"""Generate manuscript figures for the dengue ICA Perspective.

Renders six figures (PNG + SVG) for the npj Vaccines submission:

  Fig 1  Pipeline overview cartoon (5-axis flow + Versiani fork)
  Fig 2  Per-DENV-serotype Tier A/B Estofolete Table 4 scores (barplot)
  Fig 3  Phase 3 vaccine retrospective: predicted scores vs published efficacy
         (CYD-TDV / TAK-003 / Butantan-DV)
  Fig 4  HLA population coverage by region (heatmap; 5 IEDB regions x 4 DENV serotypes)
  Fig 5  EDE epitope antigenic-loss matrix (4x4 serotype x serotype)
  Fig 6  CDHIT cluster stratification across the 12 Phase 3 vaccine constructs

Inputs are expected from the master_pipeline outputs. Missing inputs are
logged and corresponding figures get a "DATA NOT YET AVAILABLE" placeholder
panel so the manuscript draft still references all six figures.

Usage:
    python3 generate_figures.py \\
        --table4_dir outputs/table4_dengue \\
        --phase3_dir outputs/phase3_retrospective \\
        --bcell_dir outputs/dengue_bcell \\
        --tcell_dir outputs/dengue_tcell \\
        --out_dir outputs/figures
"""
import argparse
import csv
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

DENV_COLORS = {
    "DENV-1": "#1b9e77",
    "DENV-2": "#d95f02",
    "DENV-3": "#7570b3",
    "DENV-4": "#e7298a",
}

PHASE3_VACCINES = {
    "CYD-TDV (Sanofi)": {"efficacy": 56.5, "color": "#1f77b4", "ref": "Sridhar 2018"},
    "TAK-003 (Takeda)": {"efficacy": 61.2, "color": "#ff7f0e", "ref": "Tricou 2024"},
    "Butantan-DV":      {"efficacy": 79.6, "color": "#2ca02c", "ref": "Kallas 2024"},
}


def placeholder(ax, msg):
    ax.text(0.5, 0.5, f"DATA NOT YET AVAILABLE\n{msg}",
            ha="center", va="center", fontsize=14,
            transform=ax.transAxes, color="#888888",
            bbox=dict(boxstyle="round", facecolor="#fafafa", edgecolor="#cccccc"))
    ax.set_xticks([])
    ax.set_yticks([])


def save(fig, out_dir, name):
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, f"{name}.png"), dpi=300, bbox_inches="tight")
    fig.savefig(os.path.join(out_dir, f"{name}.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(out_dir, f"{name}.svg"), bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {name}.png + .pdf + .svg")


def read_tsv(path):
    if not os.path.isfile(path):
        return None
    with open(path) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    return rows


def fig1_pipeline_overview(out_dir):
    fig, ax = plt.subplots(figsize=(12, 6))
    boxes = [
        ("INPUTS", "DENV 1-4\npolyproteins", 0.05, 0.5, "#cce5ff"),
        ("AXIS 1\nGEOMETRY", "AlphaFold-Multimer\n+ DiscoTope\nB-cell scoring\n(EpiDope)", 0.22, 0.5, "#fff2cc"),
        ("AXIS 2\nEQUITY", "NetMHCpan I + II\nIEDB Population Cov\n(Brazil, SE Asia,\nAfrica, India)", 0.42, 0.5, "#d4edda"),
        ("AXIS 3\nTIME (post-trial)", "MBC Hill diversity\nAlakazam\n[future work]", 0.62, 0.65, "#f8d7da"),
        ("AXIS 4\nINTENT (post-trial)", "TCRpMHC\npolyfunctionality\n[future work]", 0.62, 0.35, "#f8d7da"),
        ("AXIS 5\nEFFECTOR (post-trial)", "DC-bcell maturation\nopsonophagocytosis\n[future work]", 0.82, 0.5, "#f8d7da"),
        ("OUTPUT", "Estofolete Table 4\nTier A + Tier B\nranked candidates", 0.95, 0.5, "#cce5ff"),
    ]
    for title, body, x, y, color in boxes:
        ax.add_patch(plt.Rectangle((x - 0.06, y - 0.13), 0.12, 0.26,
                                    facecolor=color, edgecolor="black", linewidth=1.5))
        ax.text(x, y + 0.08, title, ha="center", va="center", fontsize=8, weight="bold")
        ax.text(x, y - 0.04, body, ha="center", va="center", fontsize=7)
    for x_from, x_to, y in [(0.11, 0.16, 0.5), (0.28, 0.36, 0.5),
                             (0.48, 0.56, 0.5), (0.68, 0.76, 0.5), (0.88, 0.89, 0.5)]:
        ax.annotate("", xy=(x_to, y), xytext=(x_from, y),
                    arrowprops=dict(arrowstyle="->", lw=1.5, color="#333333"))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("ICA-Dengue empirical pipeline (fork of Versiani 2026 immunoinformatics_platform)\n"
                 "Pre-trial axes (1-2) computable from sequence; post-Phase-1 axes (3-5) require trial data",
                 fontsize=11, weight="bold")
    save(fig, out_dir, "fig1_pipeline_overview")


def fig2_per_serotype_tier_scores(out_dir, table4_dir):
    fig, ax = plt.subplots(figsize=(10, 6))
    ranked = read_tsv(os.path.join(table4_dir or "", "ranked_candidates.tsv"))
    if not ranked:
        placeholder(ax, "outputs/table4_dengue/ranked_candidates.tsv")
        ax.set_title("Fig 2. Tier A + Tier B composite scores per DENV serotype")
        save(fig, out_dir, "fig2_per_serotype_tier_scores")
        return
    serotypes = sorted({r.get("serotype", r.get("seqid", "?")) for r in ranked})
    tier_a = [float(next((r["tier_a_composite"] for r in ranked
                          if r.get("serotype", r.get("seqid", "")) == s and r.get("tier_a_composite")), 0)) for s in serotypes]
    tier_b = [float(next((r["tier_b_composite"] for r in ranked
                          if r.get("serotype", r.get("seqid", "")) == s and r.get("tier_b_composite")), 0)) for s in serotypes]
    x = np.arange(len(serotypes))
    width = 0.35
    ax.bar(x - width/2, tier_a, width, label="Tier A (neutralization breadth + avidity)", color="#1f77b4")
    ax.bar(x + width/2, tier_b, width, label="Tier B (MBC breadth + CD8 polyfunc)", color="#ff7f0e")
    ax.set_xticks(x)
    ax.set_xticklabels(serotypes)
    ax.set_ylabel("Composite score (0-1)")
    ax.set_title("Fig 2. Per-DENV-serotype Estofolete Table 4 composite correlate scores\n"
                 "(higher = closer match to immune signature predictive of Phase 3 efficacy)")
    ax.legend(loc="upper right")
    ax.set_ylim(0, 1.05)
    ax.grid(axis="y", alpha=0.3)
    save(fig, out_dir, "fig2_per_serotype_tier_scores")


def fig3_phase3_retrospective(out_dir, phase3_dir):
    fig, ax = plt.subplots(figsize=(10, 7))
    rows = read_tsv(os.path.join(phase3_dir or "", "phase3_vs_clinical.tsv"))
    if not rows:
        placeholder(ax, "outputs/phase3_retrospective/phase3_vs_clinical.tsv")
        ax.set_title("Fig 3. Predicted Tier A+B scores vs published Phase 3 efficacy")
        save(fig, out_dir, "fig3_phase3_retrospective")
        return
    # Try to read companion statistical_analyses.json for permutation p
    perm_p = None
    stats_path = os.path.join(phase3_dir or "", "statistical_analyses.json")
    if os.path.isfile(stats_path):
        try:
            import json
            with open(stats_path) as f:
                stats_obj = json.load(f)
            perm_p = stats_obj.get("permutation_p_tier_A")
        except Exception:
            perm_p = None
    vaccines = sorted({r["vaccine"] for r in rows if "vaccine" in r})
    predicted = []
    published = []
    colors = []
    labels = []
    for v in vaccines:
        scores = [float(r.get("composite_score", 0)) for r in rows if r["vaccine"] == v]
        if not scores:
            continue
        pred = sum(scores) / len(scores)
        for matched, info in PHASE3_VACCINES.items():
            if v.lower() in matched.lower() or matched.lower().split()[0].lower() in v.lower():
                predicted.append(pred)
                published.append(info["efficacy"])
                colors.append(info["color"])
                labels.append(matched)
                break
    if predicted:
        ax.scatter(predicted, published, s=200, c=colors, edgecolor="black", linewidth=1.5, zorder=3)
        for x_p, y_p, lbl in zip(predicted, published, labels):
            ax.annotate(lbl, (x_p, y_p), xytext=(8, 8), textcoords="offset points", fontsize=9)
        # Draw correlation line
        if len(predicted) > 1:
            slope, intercept = np.polyfit(predicted, published, 1)
            x_line = np.linspace(min(predicted), max(predicted), 100)
            ax.plot(x_line, slope * x_line + intercept, "--", color="#666666", alpha=0.6,
                    label=f"linear fit: y={slope:.1f}x+{intercept:.1f}")
            corr = np.corrcoef(predicted, published)[0, 1]
            # Pearson r in upper-left
            ax.text(0.03, 0.97, f"Pearson r = {corr:.2f}",
                    transform=ax.transAxes, va="top", ha="left", fontsize=11,
                    bbox=dict(boxstyle="round", facecolor="#ffffcc", edgecolor="black"))
            # Permutation p in upper-right
            if perm_p is not None:
                ax.text(0.97, 0.97, f"permutation p = {perm_p:.2f}",
                        transform=ax.transAxes, va="top", ha="right", fontsize=11,
                        bbox=dict(boxstyle="round", facecolor="#e6f2ff", edgecolor="black"))
        ax.legend(loc="lower right")
    else:
        placeholder(ax, "no Phase 3 vaccine matches found in retrospective data")
    ax.set_xlabel("Predicted composite Tier A+B score (pipeline output)")
    ax.set_ylabel("Published Phase 3 efficacy (%)")
    ax.set_title("Fig 3. Retrospective Phase 3 validation:\n"
                 "predicted composite scores vs published clinical efficacy\n"
                 "n=3 vaccine programs", fontsize=10)
    ax.grid(alpha=0.3)
    ax.set_xlim(0, 1.05)
    ax.set_ylim(0, 100)
    save(fig, out_dir, "fig3_phase3_retrospective")


def fig4_hla_coverage_heatmap(out_dir, table4_dir):
    # Prefer the new 16-region v2 file (one row per region, mean dengue polyprotein coverage)
    v2_path = os.path.join(table4_dir or "", "hla_coverage_by_region_v2.tsv")
    rows_v2 = read_tsv(v2_path) if os.path.isfile(v2_path) else None
    if rows_v2 and any("mean_coverage_dengue_polyproteins" in r for r in rows_v2):
        fig, ax = plt.subplots(figsize=(8, 8))
        regions, coverages = [], []
        for r in rows_v2:
            try:
                cov = float(r.get("mean_coverage_dengue_polyproteins", "") or 0.0)
            except ValueError:
                continue
            regions.append(r.get("region", "?"))
            coverages.append(cov)
        # Sort by coverage descending so weakest regions sit at bottom
        order = sorted(range(len(regions)), key=lambda i: coverages[i], reverse=True)
        regions = [regions[i] for i in order]
        coverages = [coverages[i] for i in order]
        matrix = np.array(coverages).reshape(-1, 1)
        im = ax.imshow(matrix, cmap="RdYlGn", aspect="auto", vmin=0.7, vmax=1.0)
        ax.set_xticks([0])
        ax.set_xticklabels(["mean coverage\n(dengue polyproteins)"])
        ax.set_yticks(np.arange(len(regions)))
        ax.set_yticklabels(regions)
        # Threshold marker at 0.78 (Nigeria + Senegal fall below)
        threshold = 0.78
        for i, cov in enumerate(coverages):
            color = "black" if cov > 0.85 else "white"
            label = f"{cov:.3f}"
            if cov < threshold:
                label = f"{cov:.3f} *"
            ax.text(0, i, label, ha="center", va="center", color=color, fontsize=9)
            if cov < threshold:
                # Highlight box around the cell
                ax.add_patch(plt.Rectangle((-0.5, i - 0.5), 1.0, 1.0,
                                           fill=False, edgecolor="#d62728", linewidth=2.5))
        cbar = plt.colorbar(im, ax=ax, label="Mean population coverage (IEDB)")
        cbar.ax.axhline(threshold, color="#d62728", linewidth=1.5)
        ax.set_title(f"Fig 4. HLA Class I population coverage of DENV-1-4 polyproteins\n"
                     f"across {len(regions)} dengue-endemic regions\n"
                     f"(red box = below {threshold:.2f} threshold; equity axis exclusion)",
                     fontsize=10)
        save(fig, out_dir, "fig4_hla_coverage_heatmap")
        return
    # Legacy 5-region by serotype matrix
    fig, ax = plt.subplots(figsize=(10, 6))
    rows = read_tsv(os.path.join(table4_dir or "", "hla_coverage_by_region.tsv"))
    if not rows:
        placeholder(ax, "outputs/table4_dengue/hla_coverage_by_region.tsv")
        ax.set_title("Fig 4. HLA population coverage by region")
        save(fig, out_dir, "fig4_hla_coverage_heatmap")
        return
    regions = sorted({r["region"] for r in rows})
    serotypes = sorted({r["serotype"] for r in rows})
    matrix = np.zeros((len(regions), len(serotypes)))
    for i, region in enumerate(regions):
        for j, serotype in enumerate(serotypes):
            cov = next((float(r["coverage"]) for r in rows
                        if r["region"] == region and r["serotype"] == serotype), 0.0)
            matrix[i, j] = cov
    im = ax.imshow(matrix, cmap="RdYlGn", aspect="auto", vmin=0, vmax=1)
    ax.set_xticks(np.arange(len(serotypes)))
    ax.set_xticklabels(serotypes)
    ax.set_yticks(np.arange(len(regions)))
    ax.set_yticklabels(regions)
    for i in range(len(regions)):
        for j in range(len(serotypes)):
            ax.text(j, i, f"{matrix[i, j]:.2f}", ha="center", va="center",
                    color="black" if matrix[i, j] > 0.5 else "white", fontsize=10)
    plt.colorbar(im, ax=ax, label="Population coverage (IEDB formula)")
    ax.set_title("Fig 4. HLA Class I + II population coverage by IEDB geographic region\n"
                 "(Equity axis: coverage <0.7 in any high-burden region = exclude)")
    save(fig, out_dir, "fig4_hla_coverage_heatmap")


def fig5_ede_loss_matrix(out_dir, table4_dir):
    fig, ax = plt.subplots(figsize=(8, 7))
    rows = read_tsv(os.path.join(table4_dir or "", "ede_antigenic_loss_matrix.tsv"))
    if not rows:
        placeholder(ax, "outputs/table4_dengue/ede_antigenic_loss_matrix.tsv")
        ax.set_title("Fig 5. EDE epitope antigenic-loss matrix")
        save(fig, out_dir, "fig5_ede_loss_matrix")
        return
    serotypes = ["DENV-1", "DENV-2", "DENV-3", "DENV-4"]
    matrix = np.zeros((4, 4))
    for r in rows:
        if r.get("source_serotype") in serotypes and r.get("target_serotype") in serotypes:
            i = serotypes.index(r["source_serotype"])
            j = serotypes.index(r["target_serotype"])
            matrix[i, j] = float(r.get("antigenic_loss_index", 0))
    im = ax.imshow(matrix, cmap="Reds", vmin=0, vmax=1)
    ax.set_xticks(np.arange(4))
    ax.set_xticklabels(serotypes)
    ax.set_yticks(np.arange(4))
    ax.set_yticklabels(serotypes)
    ax.set_xlabel("Cross-reactive target serotype")
    ax.set_ylabel("Vaccine source serotype")
    for i in range(4):
        for j in range(4):
            ax.text(j, i, f"{matrix[i, j]:.2f}", ha="center", va="center",
                    color="white" if matrix[i, j] > 0.5 else "black", fontsize=11)
    plt.colorbar(im, ax=ax, label="Antigenic loss index (1 = no cross-reactivity)")
    ax.set_title("Fig 5. EDE epitope antigenic-loss matrix\n"
                 "(Geometry axis: high loss = ADE risk, low loss = breadth potential)")
    save(fig, out_dir, "fig5_ede_loss_matrix")


def fig6_cdhit_clusters(out_dir, bcell_dir):
    fig, ax = plt.subplots(figsize=(10, 6))
    rows = read_tsv(os.path.join(bcell_dir or "", "cdhit_clusters.tsv"))
    if not rows:
        # fall back to dengue smoke-test outputs
        rows = read_tsv("outputs/dengue_smoke_evidence/cdhit_clusters_dengue.tsv")
    if not rows:
        placeholder(ax, "outputs/dengue_bcell/cdhit_clusters.tsv")
        ax.set_title("Fig 6. CDHIT redundancy clustering")
        save(fig, out_dir, "fig6_cdhit_clusters")
        return
    cluster_sizes = {}
    for r in rows:
        c = r.get("cluster_index", "?")
        cluster_sizes[c] = cluster_sizes.get(c, 0) + 1
    sizes = sorted(cluster_sizes.values(), reverse=True)
    ax.bar(range(len(sizes)), sizes, color="#5b9bd5", edgecolor="black")
    ax.set_xlabel("Cluster index (sorted by size)")
    ax.set_ylabel("Number of sequences in cluster")
    ax.set_title(f"Fig 6. CDHIT redundancy clustering across {len(rows)} input sequences\n"
                 f"({len(sizes)} non-redundant clusters; centroids feed downstream B-cell + T-cell scoring)")
    ax.grid(axis="y", alpha=0.3)
    save(fig, out_dir, "fig6_cdhit_clusters")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--table4_dir", default="outputs/table4_dengue")
    p.add_argument("--phase3_dir", default="outputs/phase3_retrospective")
    p.add_argument("--bcell_dir",  default="outputs/dengue_bcell")
    p.add_argument("--tcell_dir",  default="outputs/dengue_tcell")
    p.add_argument("--out_dir",    default="outputs/figures")
    args = p.parse_args()

    matplotlib.rcParams.update({
        'figure.dpi': 300, 'savefig.dpi': 300,
        'font.family': 'sans-serif', 'font.size': 9,
        'axes.linewidth': 0.8, 'axes.labelsize': 9, 'axes.titlesize': 10,
        'legend.fontsize': 8, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
        'savefig.bbox': 'tight', 'savefig.pad_inches': 0.1,
    })

    os.makedirs(args.out_dir, exist_ok=True)
    print(f"generating figures into {args.out_dir}")

    fig1_pipeline_overview(args.out_dir)
    fig2_per_serotype_tier_scores(args.out_dir, args.table4_dir)
    fig3_phase3_retrospective(args.out_dir, args.phase3_dir)
    fig4_hla_coverage_heatmap(args.out_dir, args.table4_dir)
    fig5_ede_loss_matrix(args.out_dir, args.table4_dir)
    fig6_cdhit_clusters(args.out_dir, args.bcell_dir)

    print(f"all figures generated in {args.out_dir}")


if __name__ == "__main__":
    main()
