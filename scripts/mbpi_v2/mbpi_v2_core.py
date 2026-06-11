from __future__ import annotations

import hashlib
import math
import sys
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import optimize, stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, confusion_matrix, roc_auc_score

MBPI_DIR = Path(__file__).resolve().parents[1] / "mbpi"
if str(MBPI_DIR) not in sys.path:
    sys.path.insert(0, str(MBPI_DIR))

from mbpi_core import (  # noqa: E402
    EPS,
    PRIMARY_COMPONENTS,
    aggregate_scores,
    bh_fdr,
    compute_matrix_mbpi,
    compute_prediction_components,
    default_paths,
    derive_effect_vectors,
    ensure_out_dir,
    heterogeneity_from_matrix,
    observed_components,
    read_clr_dataset,
    score_angle,
    score_bpg,
    score_cddr,
    score_feature,
    score_loco,
    stratified_bootstrap_indices,
    summary_stats,
)


DSAS_COMPONENTS = [
    "DSAS_pooled_auroc",
    "DSAS_loco_auroc",
    "DSAS_pd_main_fdr",
    "DSAS_fixed_effect_fdr",
    "DSAS_effect_size",
]

DSAS_WEIGHTS = {
    "DSAS_pooled_auroc": 0.50,
    "DSAS_loco_auroc": 0.10,
    "DSAS_pd_main_fdr": 0.20,
    "DSAS_fixed_effect_fdr": 0.15,
    "DSAS_effect_size": 0.05,
}

PRS_COMPONENTS = [
    "S_CDDR",
    "S_angle",
    "S_BPG",
    "S_LOCO",
    "S_feature",
    "S_heterogeneity",
]

MANUAL_PRS_WEIGHTS = {
    "S_CDDR": 0.00,
    "S_angle": 0.20,
    "S_BPG": 0.00,
    "S_LOCO": 0.00,
    "S_feature": 0.00,
    "S_heterogeneity": 0.80,
}


def v2_out_dir(root: Path | str = ".") -> Path:
    return (Path(root).resolve() / "results" / "mbpi_v2").resolve()


def command_string(argv: list[str] | None = None) -> str:
    argv = sys.argv if argv is None else argv
    return "python " + " ".join(str(a).replace("\\", "/") for a in argv)


def bounded(value: float) -> float:
    if value is None or pd.isna(value):
        return np.nan
    return float(np.clip(value, 0.0, 1.0))


def auroc_signal_component(auroc: float) -> float:
    return bounded(math.sqrt(max(0.0, float(auroc) - 0.5) / 0.5))


def q_signal_component(proportion: float, saturation: float) -> float:
    return bounded(float(proportion) / saturation)


def effect_signal_component(effect_size: float, saturation: float = 0.25) -> float:
    return bounded(float(effect_size) / saturation)


def classify_dsas(dsas: float) -> str:
    if dsas < 0.33:
        return "no usable disease signal"
    if dsas < 0.66:
        return "weak/indeterminate disease signal"
    return "usable disease signal"


def classify_v2(dsas: float, prs: float) -> str:
    if dsas < 0.33:
        return "no usable disease signal"
    if dsas < 0.66:
        return "indeterminate"
    if prs > 0.66:
        return "non-portable"
    if prs >= 0.50:
        return "cohort-conditioned"
    if prs < 0.33:
        return "portable"
    return "indeterminate"


def three_class_label(regime: str) -> str:
    if regime == "null":
        return "no usable disease signal"
    if regime == "nonportable":
        return "non-portable"
    return "portable/indeterminate"


def three_class_prediction(final_class: str) -> str:
    if final_class == "no usable disease signal":
        return final_class
    if final_class in {"non-portable", "cohort-conditioned"}:
        return "non-portable"
    return "portable/indeterminate"


def pd_main_from_matrix(X: np.ndarray, y: np.ndarray, cohort: np.ndarray) -> dict[str, float]:
    from mbpi_core import build_design

    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=int)
    cohort = np.asarray(cohort, dtype=str)
    n, p = X.shape
    design = build_design(y, cohort, include_pd=True, include_cohort=True)
    coef, *_ = np.linalg.lstsq(design, X, rcond=None)
    fitted = design @ coef
    resid = X - fitted
    df_resid = n - design.shape[1]
    xtx_inv = np.linalg.pinv(design.T @ design)
    pd_var_factor = float(xtx_inv[-1, -1])
    beta_pd = coef[-1]
    sigma2 = np.sum(resid ** 2, axis=0) / max(1, df_resid)
    se = np.sqrt(np.maximum(sigma2 * pd_var_factor, EPS))
    t_stat = beta_pd / se
    p_pd = 2.0 * stats.t.sf(np.abs(t_stat), df_resid)
    q_pd = bh_fdr(p_pd)
    sig = q_pd < 0.05
    median_sig_abs = float(np.median(np.abs(beta_pd[sig]))) if np.any(sig) else 0.0
    median_top20_abs = float(np.mean(np.sort(np.abs(beta_pd))[-min(20, p):]))
    return {
        "pd_main_count": int(np.nansum(sig)),
        "pd_main_n": int(p),
        "pd_main_proportion": float(np.nansum(sig) / p),
        "median_abs_pd_effect_sig": median_sig_abs,
        "mean_abs_pd_effect_top20": median_top20_abs,
    }


def fixed_effect_from_matrix(X: np.ndarray, y: np.ndarray, cohort: np.ndarray) -> dict[str, float]:
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=int)
    cohort = np.asarray(cohort, dtype=str)
    p_values = []
    beta_fe_values = []
    for j in range(X.shape[1]):
        betas = []
        variances = []
        for level in sorted(pd.unique(cohort)):
            mask = cohort == level
            x1 = X[mask & (y == 1), j]
            x0 = X[mask & (y == 0), j]
            if len(x1) < 2 or len(x0) < 2:
                continue
            beta = float(np.mean(x1) - np.mean(x0))
            var = float(np.var(x1, ddof=1) / len(x1) + np.var(x0, ddof=1) / len(x0))
            betas.append(beta)
            variances.append(max(var, EPS))
        if len(betas) == 0:
            p_values.append(np.nan)
            beta_fe_values.append(np.nan)
            continue
        betas = np.asarray(betas)
        weights = 1.0 / np.asarray(variances)
        beta_fe = float(np.sum(weights * betas) / np.sum(weights))
        se_fe = float(np.sqrt(1.0 / np.sum(weights)))
        z = beta_fe / max(se_fe, EPS)
        p_values.append(float(2.0 * stats.norm.sf(abs(z))))
        beta_fe_values.append(beta_fe)
    q_fe = bh_fdr(p_values)
    sig = q_fe < 0.05
    return {
        "fixed_effect_count": int(np.nansum(sig)),
        "fixed_effect_n": int(X.shape[1]),
        "fixed_effect_proportion": float(np.nansum(sig) / X.shape[1]),
        "median_abs_fixed_effect_sig": float(np.nanmedian(np.abs(np.asarray(beta_fe_values)[sig]))) if np.any(sig) else 0.0,
    }


def dsas_from_parts(parts: dict[str, float]) -> float:
    return float(
        np.nansum([DSAS_WEIGHTS[k] * parts[k] for k in DSAS_COMPONENTS])
        / np.nansum([DSAS_WEIGHTS[k] for k in DSAS_COMPONENTS if pd.notna(parts.get(k))])
    )


def signal_components_from_metrics(
    pooled_auroc: float,
    loco_auroc: float,
    pd_main_proportion: float,
    fixed_effect_proportion: float,
    median_abs_pd_effect_sig: float,
) -> dict[str, float]:
    parts = {
        "DSAS_pooled_auroc": auroc_signal_component(pooled_auroc),
        "DSAS_loco_auroc": auroc_signal_component(loco_auroc),
        "DSAS_pd_main_fdr": q_signal_component(pd_main_proportion, 0.25),
        "DSAS_fixed_effect_fdr": q_signal_component(fixed_effect_proportion, 0.20),
        "DSAS_effect_size": effect_signal_component(median_abs_pd_effect_sig, 0.25),
    }
    parts["DSAS"] = dsas_from_parts(parts)
    return parts


def observed_signal_adequacy(root: Path) -> tuple[pd.DataFrame, dict[str, float]]:
    paths = default_paths(root, v2_out_dir(root))
    _, classifier_df = observed_components(paths)
    pooled = float(classifier_df["Mean_AUROC"].median())
    loco = float(classifier_df["Mean_LOCO_AUROC"].median())
    main_path = root / "results" / "fdr" / "genus_pd_main_effects_fdr.csv"
    het_path = root / "results" / "fdr" / "genus_heterogeneity_fdr.csv"
    main = pd.read_csv(main_path)
    het = pd.read_csv(het_path)
    pd_prop = float(main["sig_q05"].astype(bool).mean())
    fe_prop = float(het["sig_FE"].astype(bool).mean())
    sig_beta = main.loc[main["sig_q05"].astype(bool), "beta_PD"].astype(float)
    median_abs = float(np.median(np.abs(sig_beta))) if len(sig_beta) else 0.0
    parts = signal_components_from_metrics(pooled, loco, pd_prop, fe_prop, median_abs)
    raw = {
        "pooled_auroc_median": pooled,
        "loco_auroc_median": loco,
        "pd_main_count": int(main["sig_q05"].astype(bool).sum()),
        "pd_main_n": int(len(main)),
        "fixed_effect_count": int(het["sig_FE"].astype(bool).sum()),
        "fixed_effect_n": int(len(het)),
        "median_abs_pd_effect_sig": median_abs,
    }
    rows = []
    source = {
        "DSAS_pooled_auroc": str(paths.pooled_performance),
        "DSAS_loco_auroc": str(paths.pooled_performance),
        "DSAS_pd_main_fdr": str(main_path),
        "DSAS_fixed_effect_fdr": str(het_path),
        "DSAS_effect_size": str(main_path),
    }
    raw_text = {
        "DSAS_pooled_auroc": f"median pooled AUROC={pooled:.12g}",
        "DSAS_loco_auroc": f"median LOCO AUROC={loco:.12g}",
        "DSAS_pd_main_fdr": f"{raw['pd_main_count']}/{raw['pd_main_n']} PD-main q<0.05",
        "DSAS_fixed_effect_fdr": f"{raw['fixed_effect_count']}/{raw['fixed_effect_n']} FE q<0.05",
        "DSAS_effect_size": f"median abs beta among PD-main q<0.05={median_abs:.12g}",
    }
    for component in DSAS_COMPONENTS:
        rows.append(
            {
                "stage": "disease_signal_adequacy",
                "component": component,
                "score": parts[component],
                "weight": DSAS_WEIGHTS[component],
                "raw_value": raw_text[component],
                "source_file": source[component],
                "notes": "Higher values indicate stronger usable disease signal.",
            }
        )
    rows.append(
        {
            "stage": "disease_signal_adequacy",
            "component": "DSAS",
            "score": parts["DSAS"],
            "weight": 1.0,
            "raw_value": "weighted DSAS",
            "source_file": "; ".join(sorted(set(source.values()))),
            "notes": classify_dsas(parts["DSAS"]),
        }
    )
    return pd.DataFrame(rows), {**parts, **raw}


def observed_portability_components(root: Path) -> tuple[pd.DataFrame, dict[str, float]]:
    paths = default_paths(root, v2_out_dir(root))
    components, _ = observed_components(paths)
    primary = components[components["component"].isin(PRS_COMPONENTS)].copy()
    rows = []
    values = {}
    for row in primary.itertuples(index=False):
        values[row.component] = float(row.score)
        rows.append(
            {
                "stage": "portability_risk",
                "component": row.component,
                "score": float(row.score),
                "raw_value": row.raw_value,
                "source_file": row.source_file,
                "notes": row.notes,
            }
        )
    return pd.DataFrame(rows), values


def compute_simulation_signal_and_risk(
    X: np.ndarray,
    y: np.ndarray,
    cohort: np.ndarray,
    feature_cols: list[str],
    seed: int,
) -> dict[str, float]:
    metrics = compute_matrix_mbpi(
        X,
        y,
        cohort,
        feature_cols,
        seed=seed,
        recompute_prediction=True,
        recompute_heterogeneity=True,
        fixed_scores={"S_domain": 0.5},
    )
    pd_main = pd_main_from_matrix(X, y, cohort)
    fe = fixed_effect_from_matrix(X, y, cohort)
    dsas = signal_components_from_metrics(
        metrics["pooled_auroc"],
        metrics["loco_auroc"],
        pd_main["pd_main_proportion"],
        fe["fixed_effect_proportion"],
        pd_main["median_abs_pd_effect_sig"],
    )
    return {**metrics, **pd_main, **fe, **dsas}


def simulation_design(portable_alpha: list[float]) -> pd.DataFrame:
    rows = [
        {
            "setting": "null_permuted_within_cohort",
            "regime": "null",
            "alpha": 0.0,
            "heterogeneity_label": "none",
            "heterogeneity_level": 0.0,
            "description": "Disease labels permuted within cohort; cohort structure preserved.",
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
                "description": "Common sparse disease-effect vector injected into PD samples in every cohort.",
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
                "description": "Cohort-specific disease-effect vectors injected at fixed total magnitude.",
            }
        )
    return pd.DataFrame(rows)


def permute_within_cohort(y: np.ndarray, cohort: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    y_perm = y.copy()
    for level in sorted(pd.unique(cohort)):
        mask = cohort == level
        y_perm[mask] = rng.permutation(y_perm[mask])
    return y_perm


def centered(vec: np.ndarray) -> np.ndarray:
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
        effect = centered(alpha * common_norm * common_unit)
        X[y == 1] += effect
    else:
        h = float(setting["heterogeneity_level"])
        cohort_units = effect_info["cohort_units"]
        for level in sorted(pd.unique(cohort)):
            mask = (cohort == level) & (y == 1)
            mixed = (1.0 - h) * common_unit + h * cohort_units[level]
            norm = np.linalg.norm(mixed)
            mixed = common_unit if norm <= 0 else mixed / norm
            X[mask] += centered(alpha * common_norm * mixed)
    X = X - X.mean(axis=1, keepdims=True)
    return X, y, cohort


def run_simulations(root: Path, out_dir: Path, reps: int, seed: int) -> pd.DataFrame:
    paths = default_paths(root, out_dir)
    df, feature_cols = read_clr_dataset(paths.dataset_clr)
    design = simulation_design([0.25, 0.50, 1.00, 1.50])
    design.to_csv(out_dir / "mbpi_v2_simulation_design.csv", index=False)
    effect_info = derive_effect_vectors(df, feature_cols)
    rng = np.random.default_rng(seed)
    rows = []
    sim_id = 0
    for _, setting in design.iterrows():
        for rep in range(reps):
            sim_id += 1
            X_sim, y_sim, cohort_sim = simulate_setting(setting, df, feature_cols, effect_info, rng)
            metrics = compute_simulation_signal_and_risk(
                X_sim,
                y_sim,
                cohort_sim,
                feature_cols,
                seed + sim_id,
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
                "true_three_class": three_class_label(setting["regime"]),
            }
            row.update(metrics)
            rows.append(row)
    return pd.DataFrame(rows)


def normalize_weights(values: np.ndarray) -> np.ndarray:
    values = np.maximum(np.asarray(values, dtype=float), 0.0)
    if not np.isfinite(values).all() or values.sum() <= 0:
        values = np.ones_like(values)
    return values / values.sum()


def sigmoid(z: np.ndarray | float) -> np.ndarray | float:
    return 1.0 / (1.0 + np.exp(-np.asarray(z)))


def fit_nonnegative_logistic(X: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray, float, float]:
    def loss(params: np.ndarray) -> float:
        intercept = params[0]
        beta = params[1:]
        z = intercept + X @ beta
        return float(np.mean(np.logaddexp(0, z) - y * z))

    start = np.r_[0.0, np.repeat(0.1, X.shape[1])]
    bounds = [(None, None)] + [(0.0, None)] * X.shape[1]
    res = optimize.minimize(loss, start, method="L-BFGS-B", bounds=bounds)
    beta = res.x[1:] if res.success else np.ones(X.shape[1])
    intercept = float(res.x[0]) if res.success else 0.0
    return normalize_weights(beta), beta, intercept, float(res.fun)


def fit_elastic_net_positive(X: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    model = LogisticRegression(
        penalty="elasticnet",
        solver="saga",
        l1_ratio=0.5,
        C=1.0,
        class_weight="balanced",
        max_iter=10000,
        random_state=20260608,
    )
    model.fit(X, y)
    raw_coef = model.coef_[0]
    positive_coef = np.maximum(raw_coef, 0.0)
    return normalize_weights(positive_coef), positive_coef, float(model.intercept_[0])


def weighted_score(df: pd.DataFrame, weights: dict[str, float] | np.ndarray) -> np.ndarray:
    if isinstance(weights, dict):
        w = np.asarray([weights[c] for c in PRS_COMPONENTS], dtype=float)
    else:
        w = np.asarray(weights, dtype=float)
    w = normalize_weights(w)
    return df[PRS_COMPONENTS].astype(float).to_numpy() @ w


def strategy_scores(df: pd.DataFrame, weight_rows: pd.DataFrame) -> np.ndarray:
    mode = str(weight_rows["score_mode"].iloc[0])
    X = df[PRS_COMPONENTS].astype(float).to_numpy()
    if mode == "logistic_probability":
        coef = np.asarray([float(weight_rows.loc[weight_rows["component"] == c, "coefficient"].iloc[0]) for c in PRS_COMPONENTS])
        intercept = float(weight_rows["intercept"].iloc[0])
        return sigmoid(intercept + X @ coef)
    weights = {row.component: float(row.weight) for row in weight_rows.itertuples(index=False)}
    return weighted_score(df, weights)


def strategy_score_from_parts(parts: dict[str, float], weight_rows: pd.DataFrame) -> float:
    one = pd.DataFrame([{c: parts[c] for c in PRS_COMPONENTS}])
    return float(strategy_scores(one, weight_rows)[0])


def calibrate_prs_weights(sim: pd.DataFrame) -> tuple[pd.DataFrame, str]:
    bench = sim[sim["regime"].isin(["portable", "nonportable"])].copy()
    X = bench[PRS_COMPONENTS].astype(float).to_numpy()
    y = (bench["regime"] == "nonportable").astype(int).to_numpy()
    strategy_weights: dict[str, np.ndarray] = {}
    details: dict[str, str] = {}
    strategy_weights["equal_weight"] = normalize_weights(np.ones(len(PRS_COMPONENTS)))
    details["equal_weight"] = "Equal weights across all PRS components."
    nn_w, nn_coef, nn_intercept, nn_loss = fit_nonnegative_logistic(X, y)
    strategy_weights["nonnegative_logistic"] = nn_w
    details["nonnegative_logistic"] = f"Non-negative logistic regression coefficients normalized; log_loss={nn_loss:.6g}."
    en_w, positive_coef, en_intercept = fit_elastic_net_positive(X, y)
    strategy_weights["elastic_net_positive"] = en_w
    details["elastic_net_positive"] = "Elastic-net logistic coefficients truncated at zero; positive coefficients used as calibrated probability terms."
    coefficients = {
        "equal_weight": normalize_weights(np.ones(len(PRS_COMPONENTS))),
        "nonnegative_logistic": nn_coef,
        "elastic_net_positive": positive_coef,
        "manual_constrained": np.asarray([MANUAL_PRS_WEIGHTS[c] for c in PRS_COMPONENTS], dtype=float),
        "heterogeneity_alone": np.asarray([1.0 if c == "S_heterogeneity" else 0.0 for c in PRS_COMPONENTS]),
    }
    intercepts = {
        "equal_weight": 0.0,
        "nonnegative_logistic": nn_intercept,
        "elastic_net_positive": en_intercept,
        "manual_constrained": 0.0,
        "heterogeneity_alone": 0.0,
    }
    score_modes = {
        "equal_weight": "weighted_mean",
        "nonnegative_logistic": "logistic_probability",
        "elastic_net_positive": "logistic_probability",
        "manual_constrained": "weighted_mean",
        "heterogeneity_alone": "weighted_mean",
    }
    details["elastic_net_positive"] += "; positive_coef=" + ";".join(
        f"{c}:{v:.6g}" for c, v in zip(PRS_COMPONENTS, positive_coef)
    )
    strategy_weights["manual_constrained"] = np.asarray([MANUAL_PRS_WEIGHTS[c] for c in PRS_COMPONENTS], dtype=float)
    details["manual_constrained"] = "Manual heterogeneity-dominant non-negative score with disease-axis geometry support."
    strategy_weights["heterogeneity_alone"] = np.asarray([1.0 if c == "S_heterogeneity" else 0.0 for c in PRS_COMPONENTS])
    details["heterogeneity_alone"] = "Benchmark, not selected unless it is the best calibrated strategy."

    rows = []
    for strategy, weights in strategy_weights.items():
        weights = normalize_weights(weights)
        tmp_rows = pd.DataFrame(
            [
                {
                    "strategy": strategy,
                    "component": component,
                    "weight": float(weight),
                    "coefficient": float(coef),
                    "intercept": float(intercepts[strategy]),
                    "score_mode": score_modes[strategy],
                }
                for component, weight, coef in zip(PRS_COMPONENTS, weights, coefficients[strategy])
            ]
        )
        score = strategy_scores(bench, tmp_rows)
        auroc = roc_auc_score(y, score)
        auprc = average_precision_score(y, score)
        for component, weight, coef in zip(PRS_COMPONENTS, weights, coefficients[strategy]):
            rows.append(
                {
                    "strategy": strategy,
                    "component": component,
                    "weight": float(weight),
                    "coefficient": float(coef),
                    "intercept": float(intercepts[strategy]),
                    "score_mode": score_modes[strategy],
                    "calibration_auroc": float(auroc),
                    "calibration_auprc": float(auprc),
                    "notes": details[strategy],
                }
            )
    weight_df = pd.DataFrame(rows)
    candidate = weight_df[~weight_df["strategy"].isin(["heterogeneity_alone"])].drop_duplicates("strategy")
    selected = str(candidate.sort_values(["calibration_auroc", "calibration_auprc"], ascending=False)["strategy"].iloc[0])
    return weight_df, selected


def apply_prs_strategies(df: pd.DataFrame, weights: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for strategy, group in weights.groupby("strategy"):
        out[f"PRS_{strategy}"] = strategy_scores(out, group)
    return out


def selected_weight_dict(weights: pd.DataFrame, selected: str) -> dict[str, float]:
    group = weights[weights["strategy"] == selected]
    return {row.component: float(row.weight) for row in group.itertuples(index=False)}


def summarize_simulations(sim: pd.DataFrame, selected_strategy: str) -> pd.DataFrame:
    rows = []
    group_cols = ["regime", "setting", "alpha", "heterogeneity_label", "heterogeneity_level"]
    for keys, group in sim.groupby(group_cols, sort=False):
        row = dict(zip(group_cols, keys))
        row["n_replicates"] = len(group)
        for metric in ["DSAS", f"PRS_{selected_strategy}", "MBPI_v2_score", "pooled_auroc", "loco_auroc", "S_heterogeneity"]:
            stats_ = summary_stats(group[metric])
            row[f"{metric}_mean"] = stats_["mean"]
            row[f"{metric}_median"] = stats_["median"]
            row[f"{metric}_q2.5"] = stats_["q2.5"]
            row[f"{metric}_q97.5"] = stats_["q97.5"]
        class_counts = group["final_class"].value_counts(normalize=True).to_dict()
        row["no_signal_fraction"] = class_counts.get("no usable disease signal", 0.0)
        row["indeterminate_fraction"] = class_counts.get("indeterminate", 0.0)
        row["portable_fraction"] = class_counts.get("portable", 0.0)
        row["cohort_conditioned_fraction"] = class_counts.get("cohort-conditioned", 0.0)
        row["nonportable_fraction"] = class_counts.get("non-portable", 0.0)
        rows.append(row)
    return pd.DataFrame(rows)


def build_ablation_summary(sim: pd.DataFrame, weights: pd.DataFrame, selected_strategy: str) -> pd.DataFrame:
    rows = []
    null_vs_signal = sim[sim["regime"].isin(["null", "portable", "nonportable"])].copy()
    y_signal = (null_vs_signal["regime"] != "null").astype(int).to_numpy()
    rows.append(
        {
            "diagnostic": "DSAS: null vs disease signal",
            "task": "null_vs_disease_signal",
            "auroc": roc_auc_score(y_signal, null_vs_signal["DSAS"]),
            "auprc": average_precision_score(y_signal, null_vs_signal["DSAS"]),
            "notes": "Acceptance target AUROC >= 0.80.",
        }
    )
    bench = sim[sim["regime"].isin(["portable", "nonportable"])].copy()
    y_nonportable = (bench["regime"] == "nonportable").astype(int).to_numpy()
    for strategy in sorted(weights["strategy"].unique()):
        col = f"PRS_{strategy}"
        rows.append(
            {
                "diagnostic": f"PRS: {strategy}",
                "task": "portable_vs_nonportable",
                "auroc": roc_auc_score(y_nonportable, bench[col]),
                "auprc": average_precision_score(y_nonportable, bench[col]),
                "notes": "PRS calibration excludes null simulations.",
            }
        )
    y_full = (sim["regime"] == "nonportable").astype(int).to_numpy()
    rows.append(
        {
            "diagnostic": "MBPI-v2 gated score",
            "task": "nonportable_vs_null_or_portable",
            "auroc": roc_auc_score(y_full, sim["MBPI_v2_score"]),
            "auprc": average_precision_score(y_full, sim["MBPI_v2_score"]),
            "notes": f"Uses DSAS multiplied by PRS_{selected_strategy}.",
        }
    )
    null = sim[sim["regime"] == "null"]
    rows.append(
        {
            "diagnostic": "Null false non-portability/cohort-conditioned rate",
            "task": "null_guardrail",
            "auroc": np.nan,
            "auprc": np.nan,
            "notes": float(null["final_class"].isin(["non-portable", "cohort-conditioned"]).mean()) if len(null) else np.nan,
        }
    )
    return pd.DataFrame(rows)


def build_confusion(sim: pd.DataFrame) -> pd.DataFrame:
    labels = ["no usable disease signal", "portable/indeterminate", "non-portable"]
    y_true = sim["true_three_class"].map(three_class_prediction)
    y_pred = sim["final_class"].map(three_class_prediction)
    mat = confusion_matrix(y_true, y_pred, labels=labels)
    rows = []
    for i, true_label in enumerate(labels):
        for j, pred_label in enumerate(labels):
            rows.append({"true_class": true_label, "predicted_class": pred_label, "n": int(mat[i, j])})
    return pd.DataFrame(rows)


def bootstrap_v2(
    root: Path,
    out_dir: Path,
    selected_rows: pd.DataFrame,
    observed_fixed: dict[str, float],
    B: int,
    seed: int,
) -> pd.DataFrame:
    paths = default_paths(root, out_dir)
    df, feature_cols = read_clr_dataset(paths.dataset_clr)
    rng = np.random.default_rng(seed)
    rows = []
    for b in range(B):
        idx = stratified_bootstrap_indices(df, rng)
        boot = df.loc[idx].reset_index(drop=True)
        X = boot[feature_cols].to_numpy(float)
        y = boot["PD"].to_numpy(int)
        cohort = boot["Cohort"].to_numpy(str)
        pd_main = pd_main_from_matrix(X, y, cohort)
        fe = fixed_effect_from_matrix(X, y, cohort)
        het = heterogeneity_from_matrix(X, y, cohort)
        from mbpi_core import disease_axis_geometry, marginal_r2_cddr

        r2 = marginal_r2_cddr(X, y, cohort)
        geom = disease_axis_geometry(X, y, cohort)
        dsas_parts = signal_components_from_metrics(
            observed_fixed["pooled_auroc_median"],
            observed_fixed["loco_auroc_median"],
            pd_main["pd_main_proportion"],
            fe["fixed_effect_proportion"],
            pd_main["median_abs_pd_effect_sig"],
        )
        prs_parts = {
            "S_CDDR": score_cddr(r2["CDDR"]),
            "S_angle": score_angle(geom["Angle_Degrees"]),
            "S_BPG": observed_fixed["S_BPG"],
            "S_LOCO": observed_fixed["S_LOCO"],
            "S_feature": observed_fixed["S_feature"],
            "S_heterogeneity": het["score"],
        }
        prs = strategy_score_from_parts(prs_parts, selected_rows)
        row = {
            "replicate": b + 1,
            "DSAS": dsas_parts["DSAS"],
            "PRS": prs,
            "MBPI_v2_score": dsas_parts["DSAS"] * prs,
            "final_class": classify_v2(dsas_parts["DSAS"], prs),
            **dsas_parts,
            **pd_main,
            **fe,
            **prs_parts,
            "R2_cohort": r2["R2_cohort"],
            "R2_disease": r2["R2_disease"],
            "CDDR": r2["CDDR"],
            "angle_degrees": geom["Angle_Degrees"],
            "bootstrap_mode": "fast_genus_and_geometry",
        }
        rows.append(row)
    reps = pd.DataFrame(rows)
    reps.to_csv(out_dir / "mbpi_v2_bootstrap_replicates.csv", index=False)
    return reps


def write_observed_v2(
    root: Path,
    out_dir: Path,
    weights: pd.DataFrame,
    selected_strategy: str,
    bootstrap: pd.DataFrame,
) -> pd.DataFrame:
    signal_df, signal_raw = observed_signal_adequacy(root)
    risk_df, risk_values = observed_portability_components(root)
    selected_rows = weights[weights["strategy"] == selected_strategy]
    w = selected_weight_dict(weights, selected_strategy)
    equal_w = {c: 1.0 / len(PRS_COMPONENTS) for c in PRS_COMPONENTS}
    prs_selected = strategy_score_from_parts(risk_values, selected_rows)
    prs_equal = float(sum(risk_values[c] * equal_w[c] for c in PRS_COMPONENTS))
    dsas = float(signal_raw["DSAS"])
    final_class = classify_v2(dsas, prs_selected)
    signal_df.to_csv(out_dir / "mbpi_v2_signal_adequacy.csv", index=False)
    risk_out = risk_df.copy()
    risk_out["selected_weight"] = risk_out["component"].map(w)
    risk_out["equal_weight"] = risk_out["component"].map(equal_w)
    risk_out.to_csv(out_dir / "mbpi_v2_portability_risk.csv", index=False)
    components = pd.concat([signal_df, risk_out], ignore_index=True, sort=False)
    components.to_csv(out_dir / "mbpi_v2_components_observed.csv", index=False)

    def ci(metric: str) -> tuple[float, float]:
        return (
            float(np.quantile(bootstrap[metric], 0.025)),
            float(np.quantile(bootstrap[metric], 0.975)),
        )

    dsas_low, dsas_high = ci("DSAS")
    prs_low, prs_high = ci("PRS")
    score_low, score_high = ci("MBPI_v2_score")
    summary = pd.DataFrame(
        [
            {
                "DSAS": dsas,
                "DSAS_ci_low": dsas_low,
                "DSAS_ci_high": dsas_high,
                "DSAS_class": classify_dsas(dsas),
                "PRS": prs_selected,
                "PRS_ci_low": prs_low,
                "PRS_ci_high": prs_high,
                "PRS_strategy": selected_strategy,
                "PRS_equal_weight": prs_equal,
                "MBPI_v2_score": dsas * prs_selected,
                "MBPI_v2_score_ci_low": score_low,
                "MBPI_v2_score_ci_high": score_high,
                "final_class": final_class,
                "interpretation": "Detectable disease signal with cohort-conditioned non-portability support."
                if final_class in {"non-portable", "cohort-conditioned"}
                else final_class,
            }
        ]
    )
    summary.to_csv(out_dir / "mbpi_v2_observed_summary.csv", index=False)
    return summary


def write_algorithm_spec(root: Path, out_dir: Path, command: str, selected_strategy: str) -> None:
    spec = f"""# MBPI-v2 Algorithm Specification

MBPI-v2 is a two-stage diagnostic algorithm for microbiome biomarker portability. It separates disease-signal adequacy from portability risk. This is a technical specification, not manuscript text.

## Inputs

- CLR-transformed microbiome matrix `X` with dimension `n x p`
- phenotype vector `y`
- cohort vector `c`
- optional classifier summaries
- optional genus-level FDR and fixed-effect meta-analysis tables

Observed run primary input: `{root / 'results' / 'portability_analysis' / 'dataset_clr.csv'}`

## Stage 1: Disease Signal Adequacy

The Disease Signal Adequacy Score, DSAS, is bounded in `[0, 1]`.

Components:

- `DSAS_pooled_auroc = sqrt(max(0, median(AUROC_pooled) - 0.5) / 0.5)`
- `DSAS_loco_auroc = sqrt(max(0, median(AUROC_LOCO) - 0.5) / 0.5)`
- `DSAS_pd_main_fdr = min(1, proportion(PD-main BH q < 0.05) / 0.25)`
- `DSAS_fixed_effect_fdr = min(1, proportion(fixed-effect meta-analysis BH q < 0.05) / 0.20)`
- `DSAS_effect_size = min(1, median absolute PD beta among PD-main q < 0.05 / 0.25)`

Weights:

{pd.Series(DSAS_WEIGHTS).to_string()}

Classification:

- `DSAS < 0.33`: no usable disease signal
- `0.33 <= DSAS < 0.66`: weak/indeterminate disease signal
- `DSAS >= 0.66`: usable disease signal

## Stage 2: Conditional Portability Risk

If `DSAS < 0.33`, MBPI-v2 returns `no usable disease signal` and does not call the case non-portable.

If `0.33 <= DSAS < 0.66`, MBPI-v2 returns `indeterminate`.

If `DSAS >= 0.66`, MBPI-v2 computes PRS from:

- `S_CDDR`
- `S_angle`
- `S_BPG`
- `S_LOCO`
- `S_feature`
- `S_heterogeneity`

Equal-weight PRS and calibrated-weight PRS are both reported. For logistic calibration strategies, PRS is the calibrated logistic probability from non-negative component coefficients; normalized component weights are also reported for interpretability. For manual and equal-weight strategies, PRS is a bounded weighted mean. The primary calibrated strategy selected for this run is `{selected_strategy}`.

Final conditional rule:

- `PRS > 0.66`: non-portable
- `0.50 <= PRS <= 0.66`: cohort-conditioned
- `PRS < 0.33`: portable
- otherwise: indeterminate

## Calibration

Null simulations are retained for DSAS calibration and excluded from portable-vs-nonportable PRS calibration. PRS weights are fit using non-negative logistic regression, elastic-net logistic regression with negative coefficients truncated at zero, and a manual constrained heterogeneity-dominant score.

## Bootstrap

The bootstrap is stratified by `Cohort x PD`. The default v2 run recomputes genus-level signal proportions, fixed-effect signal proportions, CDDR, disease-axis angle, and heterogeneity from each bootstrap matrix. Existing observed classifier and feature-stability summaries are held fixed for compatibility with the observed multi-classifier analysis.

## Reproducibility

Run from the repository root:

```bash
{command}
```
"""
    (out_dir / "MBPI_v2_algorithm_specification.md").write_text(spec, encoding="utf-8")


def sha256_file(path: Path) -> str:
    if not path.exists() or not path.is_file():
        return ""
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_run_log(root: Path, out_dir: Path, command: str, selected_strategy: str, warnings: list[str]) -> None:
    observed = pd.read_csv(out_dir / "mbpi_v2_observed_summary.csv")
    ablation = pd.read_csv(out_dir / "mbpi_v2_ablation_summary.csv")
    sim_summary = pd.read_csv(out_dir / "mbpi_v2_simulation_summary.csv", keep_default_na=False)
    weights = pd.read_csv(out_dir / "mbpi_v2_weights.csv")
    paths = default_paths(root, out_dir)
    df, feature_cols = read_clr_dataset(paths.dataset_clr)
    lines = [
        "MBPI-v2 run log",
        f"date_time: {datetime.now().isoformat(timespec='seconds')}",
        f"repository_root: {root}",
        f"command: {command}",
        f"selected_prs_strategy: {selected_strategy}",
        "",
        "input_files_used:",
    ]
    for path in [
        paths.dataset_clr,
        paths.pooled_performance,
        paths.feature_stability,
        root / "results" / "fdr" / "genus_pd_main_effects_fdr.csv",
        root / "results" / "fdr" / "genus_heterogeneity_fdr.csv",
        root / "results" / "fdr" / "genus_disease_by_cohort_interactions_fdr.csv",
    ]:
        lines.append(f"- {path} | found={path.exists()} | sha256={sha256_file(path)}")
    lines.extend(
        [
            "",
            f"number_of_samples: {len(df)}",
            f"number_of_genera: {len(feature_cols)}",
            f"cohort_counts: {df.groupby('Cohort').size().to_dict()}",
            f"disease_counts: {df.groupby(['Cohort', 'PD']).size().to_dict()}",
            "",
            f"DSAS_observed: {observed['DSAS'].iloc[0]}",
            f"DSAS_95_CI: {observed['DSAS_ci_low'].iloc[0]} to {observed['DSAS_ci_high'].iloc[0]}",
            f"PRS_observed: {observed['PRS'].iloc[0]}",
            f"PRS_95_CI: {observed['PRS_ci_low'].iloc[0]} to {observed['PRS_ci_high'].iloc[0]}",
            f"MBPI_v2_final_class: {observed['final_class'].iloc[0]}",
            "",
            "selected_weights:",
        ]
    )
    sel = weights[weights["strategy"] == selected_strategy]
    for row in sel.itertuples(index=False):
        lines.append(f"- {row.component}: {row.weight}")
    lines.extend(["", "simulation_summary:"])
    keep_cols = [c for c in ["regime", "setting", "n_replicates", "DSAS_median", f"PRS_{selected_strategy}_median", "nonportable_fraction"] if c in sim_summary.columns]
    lines.extend(sim_summary[keep_cols].to_string(index=False).splitlines())
    lines.extend(["", "ablation_summary:"])
    lines.extend(ablation.to_string(index=False).splitlines())
    lines.extend(["", "warnings_or_limitations:"])
    if warnings:
        lines.extend(f"- {w}" for w in warnings)
    else:
        lines.append("- none")
    (out_dir / "mbpi_v2_run_log.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def save_pair(fig: plt.Figure, out_dir: Path, stem: str) -> None:
    fig.savefig(out_dir / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(out_dir / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def make_figures(out_dir: Path, selected_strategy: str) -> None:
    sns.set_theme(style="whitegrid")
    observed = pd.read_csv(out_dir / "mbpi_v2_observed_summary.csv")
    sim = pd.read_csv(out_dir / "mbpi_v2_simulation_results.csv", keep_default_na=False)
    boot = pd.read_csv(out_dir / "mbpi_v2_bootstrap_replicates.csv")
    ablation = pd.read_csv(out_dir / "mbpi_v2_ablation_summary.csv")

    fig, ax = plt.subplots(figsize=(10, 5.6))
    ax.set_axis_off()
    boxes = [
        ("Inputs\nCLR X, phenotype y,\ncohort c, summaries", 0.06, 0.68),
        ("Stage 1\nDisease Signal\nAdequacy Score", 0.36, 0.68),
        ("Gate\nNo signal, weak signal,\nor usable signal", 0.66, 0.68),
        ("Stage 2\nPortability Risk\nScore", 0.36, 0.28),
        ("Calibration\nSimulation weights\nand guardrails", 0.06, 0.28),
        ("Final class\nportable, indeterminate,\nnon-portable, or no signal", 0.66, 0.28),
    ]
    for text, x, y in boxes:
        ax.add_patch(plt.Rectangle((x, y), 0.24, 0.18, facecolor="#f7f7f7", edgecolor="#333333", linewidth=1.1))
        ax.text(x + 0.12, y + 0.09, text, ha="center", va="center", fontsize=10)
    for start, end in [((0.30, 0.77), (0.36, 0.77)), ((0.60, 0.77), (0.66, 0.77)), ((0.78, 0.68), (0.78, 0.46)), ((0.60, 0.37), (0.66, 0.37)), ((0.30, 0.37), (0.36, 0.37)), ((0.48, 0.68), (0.48, 0.46))]:
        ax.annotate("", xy=end, xytext=start, arrowprops=dict(arrowstyle="->", lw=1.1, color="#333333"))
    ax.set_title("MBPI-v2 Two-Stage Diagnostic Algorithm", fontsize=14, pad=10)
    save_pair(fig, out_dir, "fig_mbpi_v2_two_stage_algorithm")

    fig, ax = plt.subplots(figsize=(8.5, 6))
    palette = {"null": "#7f7f7f", "portable": "#4c78a8", "nonportable": "#c44e52"}
    sns.scatterplot(data=sim, x="DSAS", y=f"PRS_{selected_strategy}", hue="regime", palette=palette, alpha=0.55, s=35, ax=ax)
    ax.scatter(observed["DSAS"].iloc[0], observed["PRS"].iloc[0], marker="*", s=260, color="#111111", label="observed")
    ax.axvline(0.33, color="#666666", linestyle="--", linewidth=1)
    ax.axvline(0.66, color="#333333", linestyle="--", linewidth=1)
    ax.axhline(0.50, color="#666666", linestyle=":", linewidth=1)
    ax.axhline(0.66, color="#333333", linestyle=":", linewidth=1)
    ax.set_xlim(0, 1.02)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("DSAS")
    ax.set_ylabel("PRS")
    ax.set_title("Disease Signal vs Portability Risk")
    ax.legend(frameon=False)
    save_pair(fig, out_dir, "fig_mbpi_v2_signal_vs_portability_plane")

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    for ax, metric, label in zip(axes, ["DSAS", "PRS"], ["DSAS", "PRS"]):
        sns.histplot(boot[metric], bins=30, color="#4c78a8", edgecolor="white", ax=ax)
        ax.axvline(observed[metric].iloc[0], color="#111111", linewidth=2, label=f"observed {observed[metric].iloc[0]:.3f}")
        ax.axvline(np.quantile(boot[metric], 0.025), color="#c44e52", linestyle="--", linewidth=1.2)
        ax.axvline(np.quantile(boot[metric], 0.975), color="#c44e52", linestyle="--", linewidth=1.2)
        ax.set_xlabel(label)
        ax.legend(frameon=False)
    fig.suptitle("MBPI-v2 Bootstrap Distributions")
    save_pair(fig, out_dir, "fig_mbpi_v2_bootstrap_distribution")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    order = ["null", "portable", "nonportable"]
    sns.boxplot(data=sim, x="regime", y="DSAS", order=order, color="#d9e6f2", ax=axes[0])
    axes[0].axhline(0.33, color="#666666", linestyle="--", linewidth=1)
    axes[0].axhline(0.66, color="#333333", linestyle="--", linewidth=1)
    axes[0].set_title("DSAS Calibration")
    axes[0].set_xlabel("")
    axes[0].set_ylabel("DSAS")
    sns.boxplot(data=sim[sim["regime"].isin(["portable", "nonportable"])], x="regime", y=f"PRS_{selected_strategy}", order=["portable", "nonportable"], color="#f1d0d0", ax=axes[1])
    axes[1].axhline(0.66, color="#333333", linestyle="--", linewidth=1)
    axes[1].set_title("PRS Calibration")
    axes[1].set_xlabel("")
    axes[1].set_ylabel("PRS")
    save_pair(fig, out_dir, "fig_mbpi_v2_simulation_calibration")

    plot_df = ablation.copy()
    plot_df["label"] = plot_df["diagnostic"]
    plot_df["metric"] = plot_df["auroc"]
    null_guardrail = "Null false non-portability/cohort-conditioned rate"
    plot_df.loc[plot_df["diagnostic"] == null_guardrail, "metric"] = plot_df.loc[
        plot_df["diagnostic"] == null_guardrail, "notes"
    ].astype(float)
    fig, ax = plt.subplots(figsize=(10, 5.7))
    plot_df = plot_df.sort_values("metric", ascending=False)
    colors = ["#2f6f9f" if selected_strategy in x or x.startswith("MBPI-v2") else "#9ecae1" for x in plot_df["label"]]
    ax.barh(plot_df["label"], plot_df["metric"], color=colors, edgecolor="#333333", linewidth=0.5)
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("AUROC, or rate for null guardrail")
    ax.set_title("MBPI-v2 Calibration and Ablation Benchmark")
    ax.invert_yaxis()
    for y, val in enumerate(plot_df["metric"]):
        ax.text(min(float(val) + 0.015, 1.0), y, f"{float(val):.3f}", va="center", fontsize=9)
    save_pair(fig, out_dir, "fig_mbpi_v2_ablation_benchmark")
