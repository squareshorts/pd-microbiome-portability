from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from mbpi_core import (
    auprc_from_scores,
    default_paths,
    ensure_out_dir,
    optimal_accuracy_at_threshold,
    roc_auc_from_scores,
)


DIAGNOSTICS = {
    "CDDR alone": ["S_CDDR"],
    "disease-axis angle alone": ["S_angle"],
    "BPG alone": ["S_BPG"],
    "LOCO AUROC alone": ["S_LOCO"],
    "feature instability alone": ["S_feature"],
    "heterogeneity proportion alone": ["S_heterogeneity"],
    "MBPI without feature stability": ["S_CDDR", "S_angle", "S_BPG", "S_LOCO", "S_heterogeneity"],
    "MBPI without heterogeneity": ["S_CDDR", "S_angle", "S_BPG", "S_LOCO", "S_feature"],
    "MBPI without prediction metrics": ["S_CDDR", "S_angle", "S_feature", "S_heterogeneity"],
    "full MBPI": ["S_CDDR", "S_angle", "S_BPG", "S_LOCO", "S_feature", "S_heterogeneity"],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run MBPI ablation benchmark on simulation outputs.")
    parser.add_argument("--root", default=".", help="Repository root.")
    parser.add_argument("--out-dir", default="results/mbpi", help="Output directory.")
    return parser.parse_args()


def diagnostic_scores(df: pd.DataFrame, components: list[str]) -> np.ndarray:
    if len(components) == 1:
        return df[components[0]].astype(float).to_numpy()
    return df[components].astype(float).mean(axis=1).to_numpy()


def main() -> None:
    args = parse_args()
    paths = default_paths(Path(args.root), Path(args.out_dir))
    ensure_out_dir(paths.out_dir)
    sim_path = paths.out_dir / "mbpi_simulation_results.csv"
    if not sim_path.exists():
        raise FileNotFoundError(f"Missing simulation results: {sim_path}")
    sim = pd.read_csv(sim_path, keep_default_na=False)
    bench = sim[sim["regime"].isin(["portable", "nonportable"])].copy()
    bench["y_true"] = (bench["regime"] == "nonportable").astype(int)

    result_rows = []
    summary_rows = []
    for diagnostic, components in DIAGNOSTICS.items():
        scores = diagnostic_scores(bench, components)
        diagnostic_df = bench[
            [
                "simulation_id",
                "regime",
                "setting",
                "replicate",
                "alpha",
                "heterogeneity_label",
                "heterogeneity_level",
                "y_true",
            ]
        ].copy()
        diagnostic_df["diagnostic"] = diagnostic
        diagnostic_df["score"] = scores
        result_rows.append(diagnostic_df)

        valid = np.isfinite(scores)
        y_true = bench.loc[valid, "y_true"].to_numpy(int)
        valid_scores = scores[valid]
        rho, rho_p = stats.spearmanr(bench.loc[valid, "heterogeneity_level"], valid_scores)
        alpha_aurocs = []
        for alpha, alpha_group in bench.loc[valid].groupby("alpha"):
            if len(np.unique(alpha_group["y_true"])) < 2:
                continue
            alpha_scores = diagnostic_scores(alpha_group, components)
            alpha_aurocs.append(roc_auc_from_scores(alpha_group["y_true"].to_numpy(int), alpha_scores))
        summary_rows.append(
            {
                "diagnostic": diagnostic,
                "components": ";".join(components),
                "n": int(valid.sum()),
                "auroc": roc_auc_from_scores(y_true, valid_scores),
                "auprc": auprc_from_scores(y_true, valid_scores),
                "accuracy_threshold_0.66": optimal_accuracy_at_threshold(y_true, valid_scores, threshold=0.66),
                "spearman_heterogeneity": float(rho) if np.isfinite(rho) else np.nan,
                "spearman_heterogeneity_p": float(rho_p) if np.isfinite(rho_p) else np.nan,
                "min_alpha_auroc": float(np.nanmin(alpha_aurocs)) if alpha_aurocs else np.nan,
                "median_alpha_auroc": float(np.nanmedian(alpha_aurocs)) if alpha_aurocs else np.nan,
                "max_alpha_auroc": float(np.nanmax(alpha_aurocs)) if alpha_aurocs else np.nan,
            }
        )

    results = pd.concat(result_rows, ignore_index=True)
    summary = pd.DataFrame(summary_rows)
    results.to_csv(paths.out_dir / "mbpi_ablation_results.csv", index=False)
    summary.to_csv(paths.out_dir / "mbpi_ablation_summary.csv", index=False)
    print(f"Ablation results written to {paths.out_dir / 'mbpi_ablation_results.csv'}")
    print(f"Ablation summary written to {paths.out_dir / 'mbpi_ablation_summary.csv'}")


if __name__ == "__main__":
    main()
