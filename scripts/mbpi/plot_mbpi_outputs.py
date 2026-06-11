from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from mbpi_core import default_paths, ensure_out_dir


COMPONENT_LABELS = {
    "S_CDDR": "Cohort dominance",
    "S_angle": "Disease-axis misalignment",
    "S_BPG": "Prediction gap",
    "S_LOCO": "LOCO collapse",
    "S_feature": "Feature instability",
    "S_heterogeneity": "Genus heterogeneity",
    "S_domain": "Correction feasibility",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot MBPI outputs.")
    parser.add_argument("--root", default=".", help="Repository root.")
    parser.add_argument("--out-dir", default="results/mbpi", help="Output directory.")
    return parser.parse_args()


def save_pair(fig: plt.Figure, out_dir: Path, stem: str) -> None:
    fig.savefig(out_dir / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(out_dir / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_algorithm_schematic(out_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_axis_off()
    boxes = [
        ("Inputs\nCLR X, phenotype y,\ncohort c, optional summaries", 0.08, 0.72),
        ("Component diagnostics\nCDDR, angle, BPG,\nLOCO, feature stability,\nheterogeneity", 0.39, 0.72),
        ("Aggregate score\nPrimary MBPI excludes\ncorrection feasibility", 0.70, 0.72),
        ("Bootstrap uncertainty\nStratified by\ncohort x disease", 0.23, 0.34),
        ("Simulation calibration\nNull, portable,\nnon-portable regimes", 0.54, 0.34),
        ("Final output\nScore, CI, class,\naudit table", 0.70, 0.08),
    ]
    for text, x, y in boxes:
        ax.add_patch(
            plt.Rectangle(
                (x, y),
                0.23,
                0.17,
                facecolor="#f7f7f7",
                edgecolor="#333333",
                linewidth=1.2,
            )
        )
        ax.text(x + 0.115, y + 0.085, text, ha="center", va="center", fontsize=10)
    arrows = [
        ((0.31, 0.805), (0.39, 0.805)),
        ((0.62, 0.805), (0.70, 0.805)),
        ((0.505, 0.72), (0.345, 0.51)),
        ((0.505, 0.72), (0.655, 0.51)),
        ((0.655, 0.34), (0.77, 0.25)),
        ((0.345, 0.34), (0.72, 0.18)),
    ]
    for start, end in arrows:
        ax.annotate("", xy=end, xytext=start, arrowprops=dict(arrowstyle="->", lw=1.1, color="#333333"))
    ax.set_title("Microbiome Biomarker Portability Index (MBPI)", fontsize=14, pad=12)
    save_pair(fig, out_dir, "fig_mbpi_algorithm_schematic")


def plot_components(out_dir: Path) -> None:
    comp = pd.read_csv(out_dir / "mbpi_components_observed.csv")
    obs = pd.read_csv(out_dir / "mbpi_observed_summary.csv")
    primary = comp[comp["include_primary_mbpi"].astype(bool)].copy()
    primary["label"] = primary["component"].map(COMPONENT_LABELS)
    mbpi = float(obs["MBPI"].iloc[0])
    mbpi_plus = float(obs["MBPI_plus"].iloc[0])
    fig, ax = plt.subplots(figsize=(10, 5.4))
    colors = ["#4c78a8" if v <= 0.66 else "#c44e52" for v in primary["score"]]
    ax.bar(primary["label"], primary["score"], color=colors, edgecolor="#333333", linewidth=0.6)
    ax.axhline(0.33, color="#666666", linestyle="--", linewidth=1, label="0.33")
    ax.axhline(0.66, color="#333333", linestyle="--", linewidth=1, label="0.66")
    ax.axhline(mbpi, color="#111111", linewidth=2, label=f"MBPI = {mbpi:.3f}")
    ax.axhline(mbpi_plus, color="#7f3c8d", linewidth=2, linestyle="-.", label=f"MBPI_plus = {mbpi_plus:.3f}")
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Component score")
    ax.set_xlabel("")
    ax.set_title("Observed PD Validation Case")
    ax.tick_params(axis="x", rotation=30)
    ax.legend(frameon=False, loc="upper left")
    save_pair(fig, out_dir, "fig_mbpi_components_observed")


def plot_bootstrap(out_dir: Path) -> None:
    reps = pd.read_csv(out_dir / "mbpi_bootstrap_replicates.csv")
    summary = pd.read_csv(out_dir / "mbpi_summary.csv")
    mbpi_row = summary.loc[summary["metric"] == "MBPI"].iloc[0]
    median = float(mbpi_row["median"])
    low = float(mbpi_row["q2.5"])
    high = float(mbpi_row["q97.5"])
    fig, ax = plt.subplots(figsize=(8, 5))
    sns.histplot(reps["MBPI"], bins=32, color="#4c78a8", edgecolor="white", ax=ax)
    ax.axvline(median, color="#111111", linewidth=2, label=f"median {median:.3f}")
    ax.axvline(low, color="#c44e52", linewidth=1.5, linestyle="--", label=f"95% CI {low:.3f}-{high:.3f}")
    ax.axvline(high, color="#c44e52", linewidth=1.5, linestyle="--")
    ax.axvline(0.33, color="#666666", linewidth=1, linestyle=":")
    ax.axvline(0.66, color="#666666", linewidth=1, linestyle=":")
    ax.set_xlabel("MBPI")
    ax.set_ylabel("Bootstrap replicates")
    ax.set_title("Bootstrap Distribution of MBPI")
    ax.legend(frameon=False)
    save_pair(fig, out_dir, "fig_mbpi_bootstrap_distribution")


def plot_simulation(out_dir: Path) -> None:
    sim = pd.read_csv(out_dir / "mbpi_simulation_results.csv", keep_default_na=False)
    sim = sim.copy()
    sim["display"] = sim["regime"]
    sim.loc[sim["regime"] == "portable", "display"] = "portable\nalpha=" + sim.loc[sim["regime"] == "portable", "alpha"].map(lambda x: f"{x:.2f}")
    sim.loc[sim["regime"] == "nonportable", "display"] = "non-portable\n" + sim.loc[sim["regime"] == "nonportable", "heterogeneity_label"].astype(str)
    order = (
        ["null"]
        + [f"portable\nalpha={a:.2f}" for a in sorted(sim.loc[sim["regime"] == "portable", "alpha"].unique())]
        + [f"non-portable\n{x}" for x in ["low", "moderate", "high"]]
    )
    fig, ax = plt.subplots(figsize=(11, 5.5))
    sns.boxplot(data=sim, x="display", y="MBPI", order=order, color="#d9e6f2", linewidth=1, fliersize=1.5, ax=ax)
    sns.stripplot(data=sim, x="display", y="MBPI", order=order, color="#333333", size=2, alpha=0.35, ax=ax)
    ax.axhline(0.33, color="#666666", linestyle="--", linewidth=1)
    ax.axhline(0.66, color="#333333", linestyle="--", linewidth=1)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("Simulation regime")
    ax.set_ylabel("MBPI")
    ax.set_title("Simulation Calibration")
    save_pair(fig, out_dir, "fig_mbpi_simulation_calibration")


def plot_ablation(out_dir: Path) -> None:
    ablation = pd.read_csv(out_dir / "mbpi_ablation_summary.csv").sort_values("auroc", ascending=False)
    fig, ax = plt.subplots(figsize=(10, 5.8))
    colors = ["#2f6f9f" if d == "full MBPI" else "#9ecae1" for d in ablation["diagnostic"]]
    ax.barh(ablation["diagnostic"], ablation["auroc"], color=colors, edgecolor="#333333", linewidth=0.5)
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("AUROC: portable vs non-portable simulations")
    ax.set_ylabel("")
    ax.set_title("Ablation Benchmark")
    ax.invert_yaxis()
    for y, val in enumerate(ablation["auroc"]):
        ax.text(min(val + 0.015, 1.0), y, f"{val:.3f}", va="center", fontsize=9)
    save_pair(fig, out_dir, "fig_mbpi_ablation_benchmark")


def main() -> None:
    args = parse_args()
    paths = default_paths(Path(args.root), Path(args.out_dir))
    ensure_out_dir(paths.out_dir)
    sns.set_theme(style="whitegrid")
    plot_algorithm_schematic(paths.out_dir)
    plot_components(paths.out_dir)
    plot_bootstrap(paths.out_dir)
    plot_simulation(paths.out_dir)
    plot_ablation(paths.out_dir)
    print(f"MBPI figures written to {paths.out_dir}")


if __name__ == "__main__":
    main()
