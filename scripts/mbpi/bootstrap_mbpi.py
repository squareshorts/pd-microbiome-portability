from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from mbpi_core import (
    PLUS_COMPONENTS,
    PRIMARY_COMPONENTS,
    aggregate_scores,
    classify_mbpi,
    compute_matrix_mbpi,
    default_paths,
    ensure_out_dir,
    observed_components,
    stratified_bootstrap_indices,
    summary_stats,
    update_observed_summary_with_ci,
    write_observed_outputs,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Bootstrap MBPI uncertainty.")
    parser.add_argument("--root", default=".", help="Repository root.")
    parser.add_argument("--out-dir", default="results/mbpi", help="Output directory.")
    parser.add_argument("--B", type=int, default=500, help="Bootstrap replicate count.")
    parser.add_argument(
        "--mode",
        choices=["fast", "full"],
        default="fast",
        help="fast recomputes geometry and holds summary-derived components fixed; full also recomputes prediction and heterogeneity.",
    )
    parser.add_argument("--seed", type=int, default=20260608, help="Random seed.")
    return parser.parse_args()


def fixed_scores_from_observed(component_df: pd.DataFrame) -> dict[str, float]:
    return {
        row.component: float(row.score)
        for row in component_df.itertuples()
        if row.available and pd.notna(row.score)
    }


def main() -> None:
    args = parse_args()
    paths = default_paths(Path(args.root), Path(args.out_dir))
    ensure_out_dir(paths.out_dir)
    components, _ = write_observed_outputs(paths)
    observed_scores = fixed_scores_from_observed(components)
    dataset, feature_cols = pd.read_csv(paths.dataset_clr), None
    feature_cols = [c for c in dataset.columns if c not in {"SampleID", "PD", "Cohort", "Sex"}]
    dataset["PD"] = dataset["PD"].astype(int)
    dataset["Cohort"] = dataset["Cohort"].astype(str)
    X = dataset[feature_cols].to_numpy(float)

    rng = np.random.default_rng(args.seed)
    replicate_rows = []
    for b in range(args.B):
        idx = stratified_bootstrap_indices(dataset, rng)
        boot = dataset.loc[idx].reset_index(drop=True)
        Xb = boot[feature_cols].to_numpy(float)
        yb = boot["PD"].to_numpy(int)
        cb = boot["Cohort"].to_numpy(str)
        res = compute_matrix_mbpi(
            Xb,
            yb,
            cb,
            feature_cols,
            seed=args.seed + b + 1,
            recompute_prediction=args.mode == "full",
            recompute_heterogeneity=args.mode == "full",
            fixed_scores=observed_scores,
        )
        res["replicate"] = b + 1
        res["bootstrap_mode"] = args.mode
        res["n_samples"] = len(boot)
        replicate_rows.append(res)

    reps = pd.DataFrame(replicate_rows)
    reps.to_csv(paths.out_dir / "mbpi_bootstrap_replicates.csv", index=False)

    mbpi_observed = aggregate_scores(components, include_domain=False)
    plus_observed = aggregate_scores(components, include_domain=True)
    summary_rows = []
    for metric, observed_value in [("MBPI", mbpi_observed), ("MBPI_plus", plus_observed)]:
        stats = summary_stats(reps[metric])
        stats.update(
            {
                "metric": metric,
                "observed_value": observed_value,
                "bootstrap_mode": args.mode,
                "B": args.B,
                "classification": classify_mbpi(observed_value, stats["q2.5"], stats["q97.5"]),
            }
        )
        summary_rows.append(stats)
    summary = pd.DataFrame(summary_rows)[
        [
            "metric",
            "observed_value",
            "mean",
            "sd",
            "median",
            "q2.5",
            "q25",
            "q50",
            "q75",
            "q97.5",
            "n",
            "B",
            "bootstrap_mode",
            "classification",
        ]
    ]
    summary.to_csv(paths.out_dir / "mbpi_summary.csv", index=False)

    component_summary_rows = []
    for component in PLUS_COMPONENTS:
        stats = summary_stats(reps[component]) if component in reps.columns else summary_stats([])
        obs_row = components.loc[components["component"] == component]
        observed = float(obs_row["score"].iloc[0]) if not obs_row.empty and pd.notna(obs_row["score"].iloc[0]) else np.nan
        stats.update(
            {
                "component": component,
                "observed_score": observed,
                "available": bool(not obs_row.empty and obs_row["available"].iloc[0]),
                "bootstrap_mode": args.mode,
                "source_note": "recomputed" if args.mode == "full" or component in {"S_CDDR", "S_angle"} else "held fixed from observed summary",
            }
        )
        component_summary_rows.append(stats)
    component_summary = pd.DataFrame(component_summary_rows)[
        [
            "component",
            "observed_score",
            "mean",
            "sd",
            "median",
            "q2.5",
            "q25",
            "q50",
            "q75",
            "q97.5",
            "n",
            "available",
            "bootstrap_mode",
            "source_note",
        ]
    ]
    component_summary.to_csv(paths.out_dir / "mbpi_component_summary.csv", index=False)
    update_observed_summary_with_ci(paths)
    print(f"Bootstrap replicates written to {paths.out_dir / 'mbpi_bootstrap_replicates.csv'}")
    print(f"MBPI 95% CI: {summary.loc[summary['metric'] == 'MBPI', 'q2.5'].iloc[0]:.6f} to {summary.loc[summary['metric'] == 'MBPI', 'q97.5'].iloc[0]:.6f}")


if __name__ == "__main__":
    main()
