from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from mbpi_v2_core import (
    apply_prs_strategies,
    bootstrap_v2,
    build_ablation_summary,
    build_confusion,
    calibrate_prs_weights,
    classify_v2,
    command_string,
    ensure_out_dir,
    make_figures,
    observed_portability_components,
    observed_signal_adequacy,
    run_simulations,
    summarize_simulations,
    v2_out_dir,
    write_algorithm_spec,
    write_observed_v2,
    write_run_log,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the MBPI-v2 two-stage diagnostic pipeline.")
    parser.add_argument("--root", default=".", help="Repository root.")
    parser.add_argument("--simulation-reps", type=int, default=50, help="Replicates per simulation setting.")
    parser.add_argument("--bootstrap-B", type=int, default=500, help="Bootstrap replicate count.")
    parser.add_argument("--seed", type=int, default=20260608, help="Random seed.")
    return parser.parse_args()


def apply_final_classes(sim: pd.DataFrame, selected_strategy: str) -> pd.DataFrame:
    sim = sim.copy()
    sim["PRS"] = sim[f"PRS_{selected_strategy}"]
    sim["MBPI_v2_score"] = sim["DSAS"] * sim["PRS"]
    sim["final_class"] = [classify_v2(d, p) for d, p in zip(sim["DSAS"], sim["PRS"])]
    return sim


def main() -> None:
    args = parse_args()
    root = Path(args.root).resolve()
    out_dir = v2_out_dir(root)
    ensure_out_dir(out_dir)
    command = command_string()

    print("Running MBPI-v2 observed component setup", flush=True)
    signal_df, signal_raw = observed_signal_adequacy(root)
    risk_df, risk_values = observed_portability_components(root)
    signal_df.to_csv(out_dir / "mbpi_v2_signal_adequacy.csv", index=False)
    risk_df.to_csv(out_dir / "mbpi_v2_portability_risk.csv", index=False)

    print("Running MBPI-v2 simulation calibration", flush=True)
    sim = run_simulations(root, out_dir, args.simulation_reps, args.seed)
    weights, selected_strategy = calibrate_prs_weights(sim)
    weights.to_csv(out_dir / "mbpi_v2_weights.csv", index=False)
    sim = apply_prs_strategies(sim, weights)
    sim = apply_final_classes(sim, selected_strategy)
    sim.to_csv(out_dir / "mbpi_v2_simulation_results.csv", index=False)
    sim_summary = summarize_simulations(sim, selected_strategy)
    sim_summary.to_csv(out_dir / "mbpi_v2_simulation_summary.csv", index=False)

    print("Running MBPI-v2 bootstrap uncertainty", flush=True)
    selected_rows = weights[weights["strategy"] == selected_strategy]
    fixed = dict(signal_raw)
    fixed.update(risk_values)
    boot = bootstrap_v2(root, out_dir, selected_rows, fixed, args.bootstrap_B, args.seed)

    print("Writing MBPI-v2 observed summary", flush=True)
    observed_summary = write_observed_v2(root, out_dir, weights, selected_strategy, boot)

    print("Writing MBPI-v2 ablation and confusion outputs", flush=True)
    ablation = build_ablation_summary(sim, weights, selected_strategy)
    ablation.to_csv(out_dir / "mbpi_v2_ablation_summary.csv", index=False)
    confusion = build_confusion(sim)
    confusion.to_csv(out_dir / "mbpi_v2_confusion_matrix.csv", index=False)

    print("Writing MBPI-v2 figures and documentation", flush=True)
    write_algorithm_spec(root, out_dir, command, selected_strategy)
    make_figures(out_dir, selected_strategy)

    dsas_auc = float(ablation.loc[ablation["diagnostic"] == "DSAS: null vs disease signal", "auroc"].iloc[0])
    prs_auc = float(ablation.loc[ablation["diagnostic"] == f"PRS: {selected_strategy}", "auroc"].iloc[0])
    hetero_auc = float(ablation.loc[ablation["diagnostic"] == "PRS: heterogeneity_alone", "auroc"].iloc[0])
    null_false = float(ablation.loc[ablation["diagnostic"] == "Null false non-portability/cohort-conditioned rate", "notes"].iloc[0])
    warnings = [
        "Bootstrap recomputed genus-level signal, fixed-effect signal, CDDR, disease-axis angle, and heterogeneity; observed multi-classifier prediction and feature-stability summaries were held fixed.",
    ]
    if dsas_auc < 0.80:
        warnings.append(f"DSAS calibration AUROC was {dsas_auc:.3f}, below the requested 0.80 acceptance target.")
    if prs_auc <= hetero_auc:
        warnings.append(
            f"Selected PRS AUROC was {prs_auc:.3f}; heterogeneity alone remained at {hetero_auc:.3f}."
        )
    if null_false > 0.05:
        warnings.append(f"Null-regime false non-portability rate was {null_false:.3f}.")
    write_run_log(root, out_dir, command, selected_strategy, warnings)

    print(f"MBPI-v2 complete: {out_dir}", flush=True)
    print(f"Observed DSAS: {observed_summary['DSAS'].iloc[0]:.6f}", flush=True)
    print(f"Observed PRS: {observed_summary['PRS'].iloc[0]:.6f}", flush=True)
    print(f"Observed final class: {observed_summary['final_class'].iloc[0]}", flush=True)
    print(f"DSAS AUROC: {dsas_auc:.6f}; PRS AUROC: {prs_auc:.6f}; heterogeneity alone: {hetero_auc:.6f}", flush=True)


if __name__ == "__main__":
    main()
