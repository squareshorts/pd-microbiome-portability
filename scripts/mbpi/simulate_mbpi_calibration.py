from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from mbpi_core import (
    compute_matrix_mbpi,
    default_paths,
    derive_effect_vectors,
    ensure_out_dir,
    point_class,
    stratified_bootstrap_indices,
    summary_stats,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run MBPI simulation calibration regimes.")
    parser.add_argument("--root", default=".", help="Repository root.")
    parser.add_argument("--out-dir", default="results/mbpi", help="Output directory.")
    parser.add_argument("--reps", type=int, default=50, help="Replicates per simulation setting.")
    parser.add_argument("--seed", type=int, default=20260608, help="Random seed.")
    parser.add_argument(
        "--portable-alpha",
        default="0.25,0.50,1.00,1.50",
        help="Comma-separated portable effect multipliers.",
    )
    return parser.parse_args()


def simulation_design(portable_alpha: list[float]) -> pd.DataFrame:
    rows = [
        {
            "setting": "null_permuted_within_cohort",
            "regime": "null",
            "alpha": 0.0,
            "heterogeneity_label": "none",
            "heterogeneity_level": 0.0,
            "description": "Disease labels permuted within cohort; CLR matrix sampled from observed data.",
        }
    ]
    for alpha in portable_alpha:
        rows.append(
            {
                "setting": f"portable_alpha_{alpha:.2f}",
                "regime": "portable",
                "alpha": alpha,
                "heterogeneity_label": "none",
                "heterogeneity_level": 0.0,
                "description": "Common sparse cohort-adjusted disease-effect vector injected into all PD samples.",
            }
        )
    for label, level in [("low", 0.25), ("moderate", 0.50), ("high", 1.00)]:
        rows.append(
            {
                "setting": f"nonportable_{label}",
                "regime": "nonportable",
                "alpha": 1.0,
                "heterogeneity_label": label,
                "heterogeneity_level": level,
                "description": "Cohort-specific disease-effect vectors injected with fixed total magnitude.",
            }
        )
    return pd.DataFrame(rows)


def permute_within_cohort(y: np.ndarray, cohort: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    y_perm = y.copy()
    for level in sorted(pd.unique(cohort)):
        mask = cohort == level
        y_perm[mask] = rng.permutation(y_perm[mask])
    return y_perm


def centered_effect(vec: np.ndarray) -> np.ndarray:
    return vec - np.mean(vec)


def simulate_setting(
    setting: pd.Series,
    base_df: pd.DataFrame,
    feature_cols: list[str],
    effect_info: dict[str, object],
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    idx = stratified_bootstrap_indices(base_df, rng)
    sampled = base_df.loc[idx].reset_index(drop=True)
    y = sampled["PD"].to_numpy(int)
    cohort = sampled["Cohort"].to_numpy(str)
    if setting["regime"] == "null":
        X = sampled[feature_cols].to_numpy(float)
        y = permute_within_cohort(y, cohort, rng)
        return X, y, cohort

    neutral_X = effect_info["neutral_X"][idx]
    X = neutral_X.copy()
    alpha = float(setting["alpha"])
    common_unit = effect_info["common_unit"]
    common_norm = float(effect_info["common_norm"])
    if setting["regime"] == "portable":
        effect = centered_effect(alpha * common_norm * common_unit)
        X[y == 1] = X[y == 1] + effect
    else:
        h = float(setting["heterogeneity_level"])
        cohort_units = effect_info["cohort_units"]
        for level in sorted(pd.unique(cohort)):
            mask = (cohort == level) & (y == 1)
            cohort_unit = cohort_units[level]
            mixed = (1.0 - h) * common_unit + h * cohort_unit
            norm = np.linalg.norm(mixed)
            if norm <= 0:
                mixed = common_unit
                norm = np.linalg.norm(mixed)
            mixed = mixed / norm
            effect = centered_effect(alpha * common_norm * mixed)
            X[mask] = X[mask] + effect
    X = X - X.mean(axis=1, keepdims=True)
    return X, y, cohort


def summarize_simulations(results: pd.DataFrame) -> pd.DataFrame:
    rows = []
    group_cols = ["regime", "setting", "alpha", "heterogeneity_label", "heterogeneity_level"]
    for keys, group in results.groupby(group_cols, sort=False):
        row = dict(zip(group_cols, keys))
        row["n_replicates"] = len(group)
        for metric in ["MBPI", "pooled_auroc", "loco_auroc", "bpg_auroc", "S_heterogeneity", "S_feature"]:
            stats = summary_stats(group[metric])
            row[f"{metric}_mean"] = stats["mean"]
            row[f"{metric}_median"] = stats["median"]
            row[f"{metric}_q2.5"] = stats["q2.5"]
            row[f"{metric}_q97.5"] = stats["q97.5"]
        classes = group["classification"].value_counts(normalize=True).to_dict()
        row["portable_fraction"] = classes.get("portable", 0.0)
        row["indeterminate_fraction"] = classes.get("indeterminate", 0.0)
        row["nonportable_fraction"] = classes.get("non-portable", 0.0)
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    paths = default_paths(Path(args.root), Path(args.out_dir))
    ensure_out_dir(paths.out_dir)
    portable_alpha = [float(x.strip()) for x in args.portable_alpha.split(",") if x.strip()]
    design = simulation_design(portable_alpha)
    design.to_csv(paths.out_dir / "mbpi_simulation_design.csv", index=False)

    df = pd.read_csv(paths.dataset_clr)
    df["PD"] = df["PD"].astype(int)
    df["Cohort"] = df["Cohort"].astype(str)
    feature_cols = [c for c in df.columns if c not in {"SampleID", "PD", "Cohort", "Sex"}]
    effect_info = derive_effect_vectors(df, feature_cols)
    rng = np.random.default_rng(args.seed)
    rows = []
    sim_id = 0
    for setting_idx, setting in design.iterrows():
        for rep in range(args.reps):
            sim_id += 1
            X_sim, y_sim, cohort_sim = simulate_setting(setting, df, feature_cols, effect_info, rng)
            metrics = compute_matrix_mbpi(
                X_sim,
                y_sim,
                cohort_sim,
                feature_cols,
                seed=args.seed + sim_id,
                recompute_prediction=True,
                recompute_heterogeneity=True,
                fixed_scores={"S_domain": 0.5},
            )
            row = {
                "simulation_id": sim_id,
                "setting": setting["setting"],
                "regime": setting["regime"],
                "replicate": rep + 1,
                "alpha": setting["alpha"],
                "heterogeneity_label": setting["heterogeneity_label"],
                "heterogeneity_level": setting["heterogeneity_level"],
                "n_samples": len(y_sim),
                "n_features": len(feature_cols),
                "classification": point_class(metrics["MBPI"]),
            }
            row.update(metrics)
            rows.append(row)
    results = pd.DataFrame(rows)
    results.to_csv(paths.out_dir / "mbpi_simulation_results.csv", index=False)
    summary = summarize_simulations(results)
    summary.to_csv(paths.out_dir / "mbpi_simulation_summary.csv", index=False)
    print(f"Simulation design written to {paths.out_dir / 'mbpi_simulation_design.csv'}")
    print(f"Simulation results written to {paths.out_dir / 'mbpi_simulation_results.csv'}")
    print(f"Simulation summary written to {paths.out_dir / 'mbpi_simulation_summary.csv'}")


if __name__ == "__main__":
    main()
