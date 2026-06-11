from __future__ import annotations

import hashlib
import math
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    precision_recall_curve,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


PRIMARY_COMPONENTS = [
    "S_CDDR",
    "S_angle",
    "S_BPG",
    "S_LOCO",
    "S_feature",
    "S_heterogeneity",
]
PLUS_COMPONENTS = PRIMARY_COMPONENTS + ["S_domain"]
NON_FEATURE_COLUMNS = {"SampleID", "PD", "Cohort", "Sex"}
EPS = 1e-6


@dataclass
class MbpiPaths:
    root: Path
    out_dir: Path
    dataset_clr: Path
    pooled_performance: Path
    loco_performance: Path
    feature_stability: Path
    cddr: Path
    geometry: Path
    interaction_fdr: Path
    heterogeneity_fdr: Path
    domain_log: Path
    portability_gap: Path
    existing_figures: Path


def default_paths(root: Path | str = ".",
                  out_dir: Path | str = "results/mbpi") -> MbpiPaths:
    root = Path(root).resolve()
    return MbpiPaths(
        root=root,
        out_dir=(root / out_dir).resolve(),
        dataset_clr=root / "results" / "portability_analysis" / "dataset_clr.csv",
        pooled_performance=root / "results" / "portability_analysis" / "model_performance_pooled_cv.csv",
        loco_performance=root / "results" / "portability_analysis" / "model_performance_leave_one_cohort_out.csv",
        feature_stability=root / "results" / "portability_analysis" / "feature_stability_summary.csv",
        cddr=root / "results" / "portability_analysis" / "cohort_disease_dominance_ratio.csv",
        geometry=root / "results" / "portability_analysis" / "disease_vector_geometry.csv",
        interaction_fdr=root / "results" / "fdr" / "genus_disease_by_cohort_interactions_fdr.csv",
        heterogeneity_fdr=root / "results" / "fdr" / "genus_heterogeneity_fdr.csv",
        domain_log=root / "results" / "portability_analysis" / "debiasm_run_log_or_failure_reason.txt",
        portability_gap=root / "results" / "portability_analysis" / "biomarker_portability_gap.csv",
        existing_figures=root / "results" / "portability_analysis",
    )


def ensure_out_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def sha256_file(path: Path) -> str:
    if not path.exists() or not path.is_file():
        return ""
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_clr_dataset(path: Path) -> tuple[pd.DataFrame, list[str]]:
    df = pd.read_csv(path)
    required = {"SampleID", "PD", "Cohort"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")
    feature_cols = [c for c in df.columns if c not in NON_FEATURE_COLUMNS]
    if not feature_cols:
        raise ValueError(f"{path} does not contain CLR feature columns")
    df = df.copy()
    df["PD"] = df["PD"].astype(int)
    df["Cohort"] = df["Cohort"].astype(str)
    return df, feature_cols


def _cohort_dummies(cohort: Iterable[str]) -> tuple[np.ndarray, list[str]]:
    cohort = np.asarray(list(cohort), dtype=str)
    levels = sorted(pd.unique(cohort))
    dummies = []
    names = []
    for level in levels[1:]:
        dummies.append((cohort == level).astype(float))
        names.append(f"Cohort[{level}]")
    if dummies:
        return np.column_stack(dummies), names
    return np.zeros((len(cohort), 0)), names


def build_design(y: Iterable[int],
                 cohort: Iterable[str],
                 include_pd: bool = True,
                 include_cohort: bool = True,
                 include_interactions: bool = False) -> np.ndarray:
    y = np.asarray(list(y), dtype=float)
    cohort = np.asarray(list(cohort), dtype=str)
    parts = [np.ones((len(y), 1))]
    cohort_mat, _ = _cohort_dummies(cohort)
    if include_cohort and cohort_mat.shape[1]:
        parts.append(cohort_mat)
    if include_pd:
        parts.append(y.reshape(-1, 1))
    if include_interactions and cohort_mat.shape[1]:
        parts.append(cohort_mat * y.reshape(-1, 1))
    return np.column_stack(parts)


def fit_model_ss(X: np.ndarray, design: np.ndarray) -> float:
    coef, *_ = np.linalg.lstsq(design, X, rcond=None)
    fitted = design @ coef
    centered = fitted - X.mean(axis=0, keepdims=True)
    return float(np.sum(centered ** 2))


def marginal_r2_cddr(X: np.ndarray,
                     y: Iterable[int],
                     cohort: Iterable[str]) -> dict[str, float]:
    X = np.asarray(X, dtype=float)
    y = np.asarray(list(y), dtype=int)
    cohort = np.asarray(list(cohort), dtype=str)
    tss = float(np.sum((X - X.mean(axis=0, keepdims=True)) ** 2))
    if tss <= 0:
        return {"R2_cohort": np.nan, "R2_disease": np.nan, "CDDR": np.nan}
    d_pd = build_design(y, cohort, include_pd=True, include_cohort=False)
    d_cohort = build_design(y, cohort, include_pd=False, include_cohort=True)
    d_full = build_design(y, cohort, include_pd=True, include_cohort=True)
    ss_pd = fit_model_ss(X, d_pd)
    ss_cohort = fit_model_ss(X, d_cohort)
    ss_full = fit_model_ss(X, d_full)
    r2_cohort = max(0.0, (ss_full - ss_pd) / tss)
    r2_disease = max(0.0, (ss_full - ss_cohort) / tss)
    cddr = np.inf if r2_disease <= 0 else r2_cohort / r2_disease
    return {"R2_cohort": r2_cohort, "R2_disease": r2_disease, "CDDR": cddr}


def disease_axis_geometry(X: np.ndarray,
                          y: Iterable[int],
                          cohort: Iterable[str]) -> dict[str, float]:
    X = np.asarray(X, dtype=float)
    y = np.asarray(list(y), dtype=int)
    cohort = np.asarray(list(cohort), dtype=str)
    pca = PCA(n_components=1, svd_solver="full")
    pca.fit(X)
    pc1 = pca.components_[0]
    d_full = build_design(y, cohort, include_pd=True, include_cohort=True)
    coef, *_ = np.linalg.lstsq(d_full, X, rcond=None)
    disease_vector = coef[-1]
    denom = np.linalg.norm(disease_vector) * np.linalg.norm(pc1)
    if denom <= 0:
        cosine = 0.0
    else:
        cosine = abs(float(np.dot(disease_vector, pc1) / denom))
        cosine = float(np.clip(cosine, 0.0, 1.0))
    angle = float(np.degrees(np.arccos(cosine)))
    return {
        "Cosine_Similarity": cosine,
        "Angle_Degrees": angle,
        "PC1_Variance_Explained_Pct": float(pca.explained_variance_ratio_[0] * 100.0),
    }


def bounded(value: float) -> float:
    if value is None or pd.isna(value):
        return np.nan
    return float(np.clip(value, 0.0, 1.0))


def score_cddr(cddr: float) -> float:
    if cddr is None or pd.isna(cddr):
        return np.nan
    if np.isinf(cddr):
        return 1.0
    transformed = math.log10(1.0 + max(0.0, float(cddr)))
    return bounded(transformed / (1.0 + transformed))


def score_angle(angle_degrees: float) -> float:
    if angle_degrees is None or pd.isna(angle_degrees):
        return np.nan
    return bounded(float(angle_degrees) / 90.0)


def score_bpg(pooled_auroc: float, loco_auroc: float) -> float:
    bpg = float(pooled_auroc) - float(loco_auroc)
    denom = max(EPS, float(pooled_auroc) - 0.5)
    return bounded(max(0.0, bpg) / denom)


def score_loco(loco_auroc: float) -> float:
    return bounded(1.0 - min(1.0, max(0.0, float(loco_auroc) - 0.5) / 0.5))


def score_feature(jaccard_top20: float) -> float:
    return bounded(1.0 - float(jaccard_top20))


def bh_fdr(p_values: Iterable[float]) -> np.ndarray:
    p = np.asarray(list(p_values), dtype=float)
    q = np.full_like(p, np.nan, dtype=float)
    valid = np.isfinite(p)
    if not valid.any():
        return q
    pv = p[valid]
    order = np.argsort(pv)
    ranked = pv[order]
    n = len(ranked)
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.clip(adjusted, 0.0, 1.0)
    q[valid] = out
    return q


def heterogeneity_from_tables(interaction_path: Path,
                              heterogeneity_path: Path) -> dict[str, float | str]:
    interaction_count = np.nan
    interaction_n = np.nan
    q_count = np.nan
    q_n = np.nan
    sources = []
    if interaction_path.exists():
        df = pd.read_csv(interaction_path)
        interaction_n = len(df)
        if "sig_q05" in df.columns:
            interaction_count = int(pd.Series(df["sig_q05"]).astype(bool).sum())
        elif "q_interact_joint" in df.columns:
            interaction_count = int((df["q_interact_joint"].astype(float) < 0.05).sum())
        sources.append(str(interaction_path))
    if heterogeneity_path.exists():
        df = pd.read_csv(heterogeneity_path)
        q_n = len(df)
        if "sig_Q" in df.columns:
            q_count = int(pd.Series(df["sig_Q"]).astype(bool).sum())
        elif "q_Q" in df.columns:
            q_count = int((df["q_Q"].astype(float) < 0.05).sum())
        sources.append(str(heterogeneity_path))
    scores = []
    if np.isfinite(interaction_count) and interaction_n:
        scores.append(interaction_count / interaction_n)
    if np.isfinite(q_count) and q_n:
        scores.append(q_count / q_n)
    score = max(scores) if scores else np.nan
    return {
        "score": bounded(score),
        "interaction_count": interaction_count,
        "interaction_n": interaction_n,
        "cochran_q_count": q_count,
        "cochran_q_n": q_n,
        "source": "; ".join(sources),
    }


def heterogeneity_from_matrix(X: np.ndarray,
                              y: Iterable[int],
                              cohort: Iterable[str]) -> dict[str, float]:
    X = np.asarray(X, dtype=float)
    y = np.asarray(list(y), dtype=int)
    cohort = np.asarray(list(cohort), dtype=str)
    n, p = X.shape
    d_add = build_design(y, cohort, include_pd=True, include_cohort=True)
    d_full = build_design(y, cohort, include_pd=True, include_cohort=True, include_interactions=True)
    df_diff = d_full.shape[1] - d_add.shape[1]
    df_resid = n - d_full.shape[1]
    p_interaction = []
    for j in range(p):
        xj = X[:, j]
        coef_add, *_ = np.linalg.lstsq(d_add, xj, rcond=None)
        coef_full, *_ = np.linalg.lstsq(d_full, xj, rcond=None)
        rss_add = float(np.sum((xj - d_add @ coef_add) ** 2))
        rss_full = float(np.sum((xj - d_full @ coef_full) ** 2))
        if df_diff <= 0 or df_resid <= 0 or rss_full <= 0:
            p_interaction.append(np.nan)
        else:
            f_stat = max(0.0, ((rss_add - rss_full) / df_diff) / (rss_full / df_resid))
            p_interaction.append(float(stats.f.sf(f_stat, df_diff, df_resid)))
    q_interaction = bh_fdr(p_interaction)
    interaction_count = int(np.nansum(q_interaction < 0.05))

    p_q = []
    for j in range(p):
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
        if len(betas) < 2:
            p_q.append(np.nan)
            continue
        betas = np.asarray(betas)
        weights = 1.0 / np.asarray(variances)
        beta_fe = float(np.sum(weights * betas) / np.sum(weights))
        q_stat = float(np.sum(weights * (betas - beta_fe) ** 2))
        p_q.append(float(stats.chi2.sf(q_stat, len(betas) - 1)))
    q_q = bh_fdr(p_q)
    cochran_count = int(np.nansum(q_q < 0.05))
    score = max(interaction_count / p, cochran_count / p)
    return {
        "score": bounded(score),
        "interaction_count": interaction_count,
        "interaction_n": p,
        "cochran_q_count": cochran_count,
        "cochran_q_n": p,
    }


def make_logistic_model(seed: int) -> Pipeline:
    return Pipeline(
        [
            ("scale", StandardScaler()),
            (
                "model",
                LogisticRegression(
                    solver="liblinear",
                    penalty="l2",
                    C=1.0,
                    class_weight="balanced",
                    max_iter=1000,
                    random_state=seed,
                ),
            ),
        ]
    )


def _safe_auc(y_true: np.ndarray, scores: np.ndarray) -> float:
    if len(np.unique(y_true)) < 2:
        return np.nan
    return float(roc_auc_score(y_true, scores))


def _safe_auprc(y_true: np.ndarray, scores: np.ndarray) -> float:
    if len(np.unique(y_true)) < 2:
        return np.nan
    return float(average_precision_score(y_true, scores))


def compute_prediction_components(X: np.ndarray,
                                  y: Iterable[int],
                                  cohort: Iterable[str],
                                  feature_cols: list[str],
                                  seed: int = 1,
                                  n_splits: int = 5) -> dict[str, float]:
    X = np.asarray(X, dtype=float)
    y = np.asarray(list(y), dtype=int)
    cohort = np.asarray(list(cohort), dtype=str)
    min_class = min(np.bincount(y)) if len(np.unique(y)) == 2 else 0
    splits = max(2, min(n_splits, int(min_class)))
    skf = StratifiedKFold(n_splits=splits, shuffle=True, random_state=seed)
    pooled_scores = np.full(len(y), np.nan)
    pooled_pred_labels = np.full(len(y), np.nan)
    for train_idx, test_idx in skf.split(X, y):
        model = make_logistic_model(seed)
        model.fit(X[train_idx], y[train_idx])
        scores = model.predict_proba(X[test_idx])[:, 1]
        pooled_scores[test_idx] = scores
        pooled_pred_labels[test_idx] = (scores >= 0.5).astype(int)
    pooled_auroc = _safe_auc(y, pooled_scores)
    pooled_auprc = _safe_auprc(y, pooled_scores)
    pooled_balacc = float(balanced_accuracy_score(y, pooled_pred_labels))

    loco_aurocs = []
    loco_auprcs = []
    top20_sets = []
    for level in sorted(pd.unique(cohort)):
        test = cohort == level
        train = ~test
        model = make_logistic_model(seed)
        model.fit(X[train], y[train])
        scores = model.predict_proba(X[test])[:, 1]
        loco_aurocs.append(_safe_auc(y[test], scores))
        loco_auprcs.append(_safe_auprc(y[test], scores))
        coefs = model.named_steps["model"].coef_[0] / model.named_steps["scale"].scale_
        k = min(20, len(feature_cols))
        top_idx = np.argsort(np.abs(coefs))[::-1][:k]
        top20_sets.append(set(feature_cols[i] for i in top_idx))
    loco_auroc = float(np.nanmean(loco_aurocs))
    loco_auprc = float(np.nanmean(loco_auprcs))
    pairwise_j = []
    for i in range(len(top20_sets)):
        for j in range(i + 1, len(top20_sets)):
            union = top20_sets[i] | top20_sets[j]
            pairwise_j.append(len(top20_sets[i] & top20_sets[j]) / len(union) if union else np.nan)
    jaccard = float(np.nanmean(pairwise_j)) if pairwise_j else np.nan
    bpg = float(pooled_auroc - loco_auroc)
    return {
        "pooled_auroc": pooled_auroc,
        "pooled_auprc": pooled_auprc,
        "pooled_balanced_accuracy": pooled_balacc,
        "loco_auroc": loco_auroc,
        "loco_auprc": loco_auprc,
        "bpg_auroc": bpg,
        "feature_jaccard_top20": jaccard,
        "S_BPG": score_bpg(pooled_auroc, loco_auroc),
        "S_LOCO": score_loco(loco_auroc),
        "S_feature": score_feature(jaccard),
    }


def aggregate_scores(component_rows: pd.DataFrame,
                     include_domain: bool = False) -> float:
    names = PLUS_COMPONENTS if include_domain else PRIMARY_COMPONENTS
    rows = component_rows[
        component_rows["component"].isin(names)
        & component_rows["available"].astype(bool)
    ]
    if rows.empty:
        return np.nan
    return float(rows["score"].astype(float).mean())


def conservative_flag(component_rows: pd.DataFrame) -> bool:
    rows = component_rows[
        component_rows["component"].isin(PRIMARY_COMPONENTS)
        & component_rows["available"].astype(bool)
    ]
    return int((rows["score"].astype(float) > 0.66).sum()) >= 4


def classify_mbpi(score: float,
                  ci_low: float | None = None,
                  ci_high: float | None = None) -> str:
    if score is None or pd.isna(score):
        return "indeterminate"
    if ci_low is None or ci_high is None or pd.isna(ci_low) or pd.isna(ci_high):
        if score < 0.33:
            return "portable_point_only"
        if score > 0.66:
            return "non-portable_point_only"
        return "indeterminate_point_only"
    if score < 0.33 and ci_high < 0.50:
        return "portable"
    if score > 0.66 and ci_low > 0.50:
        return "non-portable"
    return "indeterminate"


def point_class(score: float) -> str:
    if score < 0.33:
        return "portable"
    if score > 0.66:
        return "non-portable"
    return "indeterminate"


def observed_components(paths: MbpiPaths) -> tuple[pd.DataFrame, pd.DataFrame]:
    df, feature_cols = read_clr_dataset(paths.dataset_clr)
    X = df[feature_cols].to_numpy(float)
    y = df["PD"].to_numpy(int)
    cohort = df["Cohort"].to_numpy(str)
    rows = []
    classifier_rows = []

    if paths.cddr.exists():
        cddr_df = pd.read_csv(paths.cddr)
        c_row = cddr_df.loc[cddr_df["Term"].astype(str).str.lower() == "cohort"].iloc[0]
        d_row = cddr_df.loc[cddr_df["Term"].astype(str).str.lower() == "pd"].iloc[0]
        cddr = float(c_row["CDDR"]) if "CDDR" in c_row and pd.notna(c_row["CDDR"]) else float(c_row["R2"]) / float(d_row["R2"])
        raw = f"CDDR={cddr:.12g}; R2_cohort={float(c_row['R2']):.12g}; R2_disease={float(d_row['R2']):.12g}"
        source = str(paths.cddr)
        notes = "Observed value read from existing portability output."
    else:
        r2 = marginal_r2_cddr(X, y, cohort)
        cddr = r2["CDDR"]
        raw = f"CDDR={cddr:.12g}; R2_cohort={r2['R2_cohort']:.12g}; R2_disease={r2['R2_disease']:.12g}"
        source = str(paths.dataset_clr)
        notes = "Existing CDDR summary missing; recomputed from CLR matrix."
    rows.append(_component_row("S_CDDR", score_cddr(cddr), True, raw, source, notes, True, True))

    if paths.geometry.exists():
        geom_df = pd.read_csv(paths.geometry)
        geom = dict(zip(geom_df["Metric"], geom_df["Value"]))
        angle = float(geom["Angle_Degrees"])
        raw = "; ".join(f"{k}={float(v):.12g}" for k, v in geom.items())
        source = str(paths.geometry)
        notes = "Observed disease-axis geometry read from existing output."
    else:
        geom = disease_axis_geometry(X, y, cohort)
        angle = geom["Angle_Degrees"]
        raw = "; ".join(f"{k}={float(v):.12g}" for k, v in geom.items())
        source = str(paths.dataset_clr)
        notes = "Existing disease-axis geometry missing; recomputed from CLR matrix."
    rows.append(_component_row("S_angle", score_angle(angle), True, raw, source, notes, True, True))

    if paths.pooled_performance.exists():
        perf = pd.read_csv(paths.pooled_performance)
        for _, r in perf.iterrows():
            model = str(r["Model"])
            pooled = float(r["Mean_AUROC"])
            loco = float(r["Mean_LOCO_AUROC"]) if "Mean_LOCO_AUROC" in perf.columns else np.nan
            bpg = float(r["BPG_AUROC"]) if "BPG_AUROC" in perf.columns else pooled - loco
            classifier_rows.append(
                {
                    "Model": model,
                    "Mean_AUROC": pooled,
                    "Mean_LOCO_AUROC": loco,
                    "BPG_AUROC": bpg,
                    "S_BPG": score_bpg(pooled, loco),
                    "S_LOCO": score_loco(loco),
                    "source": str(paths.pooled_performance),
                }
            )
        cdf = pd.DataFrame(classifier_rows)
        rows.append(
            _component_row(
                "S_BPG",
                float(cdf["S_BPG"].median()),
                True,
                "median over classifiers; " + "; ".join(f"{r.Model}={r.S_BPG:.4f}" for r in cdf.itertuples()),
                str(paths.pooled_performance),
                "Predictive portability-gap component computed from existing pooled and LOCO AUROC summary.",
                True,
                True,
            )
        )
        rows.append(
            _component_row(
                "S_LOCO",
                float(cdf["S_LOCO"].median()),
                True,
                "median over classifiers; " + "; ".join(f"{r.Model}={r.S_LOCO:.4f}" for r in cdf.itertuples()),
                str(paths.pooled_performance),
                "LOCO collapse component computed from existing mean LOCO AUROC summary.",
                True,
                True,
            )
        )
    else:
        rows.append(_component_row("S_BPG", np.nan, False, "", "", "Pooled performance summary missing.", True, True))
        rows.append(_component_row("S_LOCO", np.nan, False, "", "", "LOCO performance summary missing.", True, True))

    if paths.feature_stability.exists():
        fs = pd.read_csv(paths.feature_stability)
        if classifier_rows:
            cdf = pd.DataFrame(classifier_rows)
            cdf = cdf.merge(fs[["Model", "Mean_Jaccard_Top20"]], on="Model", how="left")
            cdf["S_feature"] = cdf["Mean_Jaccard_Top20"].map(score_feature)
            classifier_rows = cdf.to_dict("records")
        else:
            fs = fs.copy()
            fs["S_feature"] = fs["Mean_Jaccard_Top20"].map(score_feature)
        source_df = pd.DataFrame(classifier_rows) if classifier_rows else fs
        rows.append(
            _component_row(
                "S_feature",
                float(source_df["S_feature"].median()),
                True,
                "median over classifiers; " + "; ".join(
                    f"{r.Model}={r.S_feature:.4f}" for r in source_df.itertuples()
                ),
                str(paths.feature_stability),
                "Feature-instability component computed as 1 - top-20 Jaccard similarity.",
                True,
                True,
            )
        )
    else:
        rows.append(_component_row("S_feature", np.nan, False, "", "", "Feature stability summary missing.", True, True))

    het = heterogeneity_from_tables(paths.interaction_fdr, paths.heterogeneity_fdr)
    if pd.notna(het["score"]):
        raw = (
            f"PD_by_cohort={het['interaction_count']}/{het['interaction_n']}; "
            f"Cochran_Q={het['cochran_q_count']}/{het['cochran_q_n']}"
        )
        notes = "Heterogeneity component is the maximum of interaction and Cochran-Q FDR proportions."
        rows.append(_component_row("S_heterogeneity", het["score"], True, raw, str(het["source"]), notes, True, True))
    else:
        rows.append(_component_row("S_heterogeneity", np.nan, False, "", "", "Genus-level heterogeneity summaries missing.", True, True))

    domain_score = 0.5
    domain_notes = "No correction procedure was evaluable."
    domain_raw = "S_domain=0.5"
    domain_source = ""
    if paths.domain_log.exists():
        text = paths.domain_log.read_text(encoding="utf-8", errors="replace")
        domain_source = str(paths.domain_log)
        lower = text.lower()
        if "leaky/invalid" in lower or "not valid for strict loco" in lower:
            domain_score = 1.0
            domain_notes = (
                "Correction procedure was flagged as unusable under strict held-out-cohort validation "
                "without held-out phenotype-label leakage."
            )
        elif "valid" in lower and ("improves" in lower or "preserves" in lower):
            domain_score = 0.0
            domain_notes = "At least one correction procedure was valid and improved or preserved LOCO performance."
        domain_raw = "S_domain=" + str(domain_score)
    rows.append(_component_row("S_domain", domain_score, True, domain_raw, domain_source, domain_notes, False, True))

    component_df = pd.DataFrame(rows)
    classifier_df = pd.DataFrame(classifier_rows)
    return component_df, classifier_df


def _component_row(component: str,
                   score: float,
                   available: bool,
                   raw_value: str,
                   source: str,
                   notes: str,
                   include_primary: bool,
                   include_plus: bool) -> dict[str, object]:
    return {
        "component": component,
        "score": score,
        "available": bool(available),
        "raw_value": raw_value,
        "source_file": source,
        "include_primary_mbpi": include_primary,
        "include_mbpi_plus": include_plus,
        "exceeds_0.66": bool(pd.notna(score) and score > 0.66),
        "notes": notes,
    }


def validation_rows(paths: MbpiPaths) -> pd.DataFrame:
    items = [
        ("CLR matrix", paths.dataset_clr, "Primary 683 x 62 CLR feature source with SampleID, PD, and Cohort columns."),
        ("metadata", paths.dataset_clr, "Metadata columns are embedded in the primary CLR dataset."),
        ("pooled model performance", paths.pooled_performance, "Existing pooled CV and mean LOCO model summary."),
        ("LOCO model performance", paths.loco_performance, "Existing held-out-cohort performance by classifier and cohort."),
        ("portability gap summary", paths.portability_gap, "Existing BPG summary."),
        ("feature-stability summary", paths.feature_stability, "Existing top-feature overlap/Jaccard summary."),
        ("genus PD-by-cohort FDR", paths.interaction_fdr, "Existing genus-level interaction FDR table."),
        ("genus Cochran-Q FDR", paths.heterogeneity_fdr, "Existing genus-level heterogeneity FDR table."),
        ("domain-correction feasibility", paths.domain_log, "Existing DEBIAS-M feasibility log."),
        ("existing portability figures", paths.existing_figures, "Directory inspected for existing figures; not used as MBPI numeric input."),
    ]
    rows = []
    for name, path, notes in items:
        found = path.exists()
        rows_count = ""
        cols_count = ""
        file_used = str(path) if found else ""
        if found and path.is_file() and path.suffix.lower() == ".csv":
            try:
                df = pd.read_csv(path)
                rows_count = len(df)
                cols_count = len(df.columns)
            except Exception as exc:
                notes = f"{notes} Could not read shape: {exc}"
        elif found and path.is_dir():
            figure_count = len(list(path.glob("Figure_*"))) + len(list(path.glob("fig_*")))
            rows_count = figure_count
            cols_count = ""
            notes = f"{notes} Figure-like files found: {figure_count}."
        rows.append(
            {
                "required_input": name,
                "file_used": file_used,
                "found": bool(found),
                "rows": rows_count,
                "columns": cols_count,
                "notes": notes,
            }
        )
    return pd.DataFrame(rows)


def summary_stats(values: Iterable[float]) -> dict[str, float]:
    arr = np.asarray(list(values), dtype=float)
    arr = arr[np.isfinite(arr)]
    if len(arr) == 0:
        return {
            "n": 0,
            "mean": np.nan,
            "sd": np.nan,
            "median": np.nan,
            "q2.5": np.nan,
            "q25": np.nan,
            "q50": np.nan,
            "q75": np.nan,
            "q97.5": np.nan,
        }
    return {
        "n": int(len(arr)),
        "mean": float(np.mean(arr)),
        "sd": float(np.std(arr, ddof=1)) if len(arr) > 1 else 0.0,
        "median": float(np.median(arr)),
        "q2.5": float(np.quantile(arr, 0.025)),
        "q25": float(np.quantile(arr, 0.25)),
        "q50": float(np.quantile(arr, 0.50)),
        "q75": float(np.quantile(arr, 0.75)),
        "q97.5": float(np.quantile(arr, 0.975)),
    }


def write_algorithm_spec(paths: MbpiPaths, command: str) -> None:
    ensure_out_dir(paths.out_dir)
    spec = f"""# Microbiome Biomarker Portability Index (MBPI)

This technical note defines the Microbiome Biomarker Portability Index (MBPI) as a reproducible diagnostic algorithm for microbiome biomarker non-portability. It is an algorithm specification, not manuscript prose.

## Inputs

- `X`: CLR-transformed microbiome matrix with dimension `n x p`.
- `y`: binary phenotype vector, where `1` denotes PD and `0` denotes control.
- `c`: cohort vector.
- Optional classifier outputs for pooled and leave-one-cohort-out (LOCO) AUROC.
- Optional feature-importance stability outputs.
- Optional genus-level p-value/FDR tables.

The observed PD validation run uses `{paths.dataset_clr}` as the primary CLR and metadata source.

## Component Definitions

Let `R2_cohort` and `R2_disease` be marginal Aitchison-space variance fractions from the full model `X ~ Cohort + PD`, with each term evaluated conditional on the other term.

`CDDR = R2_cohort / R2_disease`

`S_CDDR = log10(1 + CDDR) / (1 + log10(1 + CDDR))`

Let `theta` be the acute angle, in degrees, between the cohort-adjusted disease-effect vector from `X ~ Cohort + PD` and PC1 of `X`.

`S_angle = theta / 90`

For each classifier `m`:

`BPG_m = AUROC_pooled,m - AUROC_LOCO,m`

`S_BPG,m = max(0, BPG_m) / max(1e-6, AUROC_pooled,m - 0.5)`

The implementation clips component scores to `[0, 1]` before aggregation so the index remains bounded. `S_BPG` is the median across classifiers.

For each classifier `m`:

`S_LOCO,m = 1 - min(1, max(0, AUROC_LOCO,m - 0.5) / 0.5)`

`S_LOCO` is the median across classifiers.

For top-20 feature Jaccard similarity `J_m`:

`S_feature,m = 1 - J_m`

`S_feature` is the median across classifiers.

For genus-level heterogeneity:

`S_heterogeneity = max(N_interaction_q05 / N_interaction_tested, N_CochranQ_q05 / N_CochranQ_tested)`

For domain-correction feasibility:

- `S_domain = 1` when the evaluated correction procedure cannot be used under strict held-out-cohort validation without held-out phenotype-label leakage.
- `S_domain = 0.5` when no correction procedure was evaluable.
- `S_domain = 0` when at least one correction procedure is valid under strict held-out-cohort validation and improves or preserves LOCO performance.

`S_domain` is reported separately and included only in `MBPI_plus`.

## Aggregation Rule

Primary MBPI excludes `S_domain`:

`MBPI = mean(S_CDDR, S_angle, S_BPG, S_LOCO, S_feature, S_heterogeneity)`

`MBPI_plus = mean(S_CDDR, S_angle, S_BPG, S_LOCO, S_feature, S_heterogeneity, S_domain)`

If a component is unavailable, the score is averaged over available components and the missing component remains in the audit table.

## Classification Rule

These are operational diagnostic thresholds for MBPI calibration, not universal biological constants.

- Portable: `MBPI < 0.33` and the 95% bootstrap CI upper bound is `< 0.50`.
- Indeterminate: `MBPI` is between `0.33` and `0.66`, or the confidence interval crosses `0.50`.
- Non-portable: `MBPI > 0.66` and the 95% bootstrap CI lower bound is `> 0.50`.

Conservative binary flag:

`non_portable_flag = TRUE` when at least four of the six primary components exceed `0.66`.

## Bootstrap Procedure

Participant-level bootstrap resampling is stratified by `Cohort x PD`. For each replicate the implementation resamples participants within each stratum, recomputes CDDR and disease-axis angle from the resampled CLR matrix, and recomputes MBPI. In fast mode, prediction, feature-stability, and genus-level heterogeneity components are held at the observed summary values; full mode can recompute those components but is slower.

Outputs:

- `mbpi_bootstrap_replicates.csv`
- `mbpi_summary.csv`
- `mbpi_component_summary.csv`

## Simulation Calibration Procedure

Simulation calibration uses the observed CLR matrix as the base and preserves cohort structure.

- Null regime: disease labels are permuted within cohort.
- Portable regime: a common sparse cohort-adjusted disease-effect vector is injected into PD samples in every cohort across an effect-size grid.
- Non-portable regime: cohort-specific disease-effect vectors are injected with low, moderate, or high heterogeneity while preserving the total effect magnitude.

Outputs:

- `mbpi_simulation_design.csv`
- `mbpi_simulation_results.csv`
- `mbpi_simulation_summary.csv`

## Ablation Procedure

The ablation benchmark compares full MBPI with individual diagnostics and reduced MBPI variants using portable vs non-portable simulation labels. Reported metrics include AUROC, AUPRC, accuracy at threshold `0.66`, Spearman correlation with simulated heterogeneity level, and robustness summaries by effect magnitude.

Outputs:

- `mbpi_ablation_results.csv`
- `mbpi_ablation_summary.csv`

## Output Field Interpretation

- `score`: bounded diagnostic score where larger values indicate stronger evidence for non-portability.
- `available`: whether the component was available for aggregation.
- `source_file`: file used to compute the observed component.
- `exceeds_0.66`: whether the component crosses the conservative component threshold.
- `MBPI`: primary score excluding domain correction feasibility.
- `MBPI_plus`: secondary score including domain correction feasibility.

## Reproducibility

Run from the repository root:

```bash
{command}
```

The run log records input files, sample counts, component values, bootstrap intervals, simulation settings, ablation summaries, warnings, and file hashes where feasible.
"""
    (paths.out_dir / "MBPI_algorithm_specification.md").write_text(spec, encoding="utf-8")


def write_run_log(paths: MbpiPaths, command: str, warnings: list[str] | None = None) -> None:
    warnings = warnings or []
    df, feature_cols = read_clr_dataset(paths.dataset_clr)
    summary_path = paths.out_dir / "mbpi_summary.csv"
    observed_path = paths.out_dir / "mbpi_observed_summary.csv"
    comp_path = paths.out_dir / "mbpi_components_observed.csv"
    sim_summary_path = paths.out_dir / "mbpi_simulation_summary.csv"
    ablation_path = paths.out_dir / "mbpi_ablation_summary.csv"
    summary = pd.read_csv(summary_path) if summary_path.exists() else pd.DataFrame()
    observed = pd.read_csv(observed_path) if observed_path.exists() else pd.DataFrame()
    comp = pd.read_csv(comp_path) if comp_path.exists() else pd.DataFrame()
    sim_summary = pd.read_csv(sim_summary_path, keep_default_na=False) if sim_summary_path.exists() else pd.DataFrame()
    ablation = pd.read_csv(ablation_path, keep_default_na=False) if ablation_path.exists() else pd.DataFrame()

    def metric_row(metric: str, col: str, default: str = "NA") -> str:
        if summary.empty or metric not in set(summary["metric"]):
            return default
        val = summary.loc[summary["metric"] == metric, col].iloc[0]
        return str(val)

    cohort_counts = df.groupby("Cohort").size().to_dict()
    disease_counts = df.groupby(["Cohort", "PD"]).size().reset_index(name="n")
    input_files = [
        paths.dataset_clr,
        paths.pooled_performance,
        paths.loco_performance,
        paths.feature_stability,
        paths.cddr,
        paths.geometry,
        paths.interaction_fdr,
        paths.heterogeneity_fdr,
        paths.domain_log,
    ]
    lines = []
    lines.append("MBPI run log")
    lines.append(f"date_time: {datetime.now().isoformat(timespec='seconds')}")
    lines.append(f"repository_root: {paths.root}")
    lines.append(f"command: {command}")
    lines.append("")
    lines.append("input_files_used:")
    for path in input_files:
        lines.append(f"- {path} | found={path.exists()} | sha256={sha256_file(path)}")
    lines.append("")
    lines.append(f"number_of_samples: {len(df)}")
    lines.append(f"number_of_genera: {len(feature_cols)}")
    lines.append(f"cohort_counts: {cohort_counts}")
    lines.append("disease_counts_by_cohort:")
    for row in disease_counts.itertuples(index=False):
        lines.append(f"- Cohort={row.Cohort}, PD={row.PD}, n={row.n}")
    lines.append("")
    lines.append(f"MBPI_value: {metric_row('MBPI', 'observed_value')}")
    lines.append(f"MBPI_95_CI: {metric_row('MBPI', 'q2.5')} to {metric_row('MBPI', 'q97.5')}")
    lines.append(f"MBPI_classification: {metric_row('MBPI', 'classification')}")
    lines.append(f"MBPI_plus_value: {metric_row('MBPI_plus', 'observed_value')}")
    lines.append(f"MBPI_plus_95_CI: {metric_row('MBPI_plus', 'q2.5')} to {metric_row('MBPI_plus', 'q97.5')}")
    lines.append(f"MBPI_plus_classification: {metric_row('MBPI_plus', 'classification')}")
    if not observed.empty and "non_portable_flag" in observed.columns:
        lines.append(f"conservative_non_portable_flag: {observed['non_portable_flag'].iloc[0]}")
    lines.append("")
    lines.append("component_scores:")
    if not comp.empty:
        for row in comp.itertuples(index=False):
            lines.append(f"- {row.component}: score={row.score}, available={row.available}, source={row.source_file}")
    lines.append("")
    lines.append("simulation_settings:")
    if not sim_summary.empty:
        settings_cols = [c for c in ["regime", "setting", "alpha", "heterogeneity_label", "n_replicates", "MBPI_mean", "MBPI_median"] if c in sim_summary.columns]
        lines.extend(sim_summary[settings_cols].to_string(index=False).splitlines())
    else:
        lines.append("NA")
    lines.append("")
    lines.append("ablation_results_summary:")
    if not ablation.empty:
        cols = [c for c in ["diagnostic", "auroc", "auprc", "accuracy_threshold_0.66", "spearman_heterogeneity"] if c in ablation.columns]
        lines.extend(ablation[cols].to_string(index=False).splitlines())
    else:
        lines.append("NA")
    lines.append("")
    lines.append("warnings_or_missing_inputs:")
    validation = validation_rows(paths)
    missing = validation.loc[~validation["found"].astype(bool), "required_input"].tolist()
    for item in missing:
        lines.append(f"- missing input: {item}")
    for warning in warnings:
        lines.append(f"- {warning}")
    if not missing and not warnings:
        lines.append("- none")
    (paths.out_dir / "mbpi_run_log.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_observed_outputs(paths: MbpiPaths) -> tuple[pd.DataFrame, pd.DataFrame]:
    ensure_out_dir(paths.out_dir)
    components, classifier_components = observed_components(paths)
    components.to_csv(paths.out_dir / "mbpi_components_observed.csv", index=False)
    components.to_csv(paths.out_dir / "mbpi_audit_table.csv", index=False)
    if not classifier_components.empty:
        classifier_components.to_csv(paths.out_dir / "mbpi_classifier_components.csv", index=False)
    validation_rows(paths).to_csv(paths.out_dir / "mbpi_input_validation.csv", index=False)
    mbpi = aggregate_scores(components, include_domain=False)
    mbpi_plus = aggregate_scores(components, include_domain=True)
    observed = pd.DataFrame(
        [
            {
                "MBPI": mbpi,
                "MBPI_plus": mbpi_plus,
                "classification": classify_mbpi(mbpi),
                "classification_plus": classify_mbpi(mbpi_plus),
                "non_portable_flag": conservative_flag(components),
                "n_primary_components_available": int(
                    components[
                        components["component"].isin(PRIMARY_COMPONENTS)
                        & components["available"].astype(bool)
                    ].shape[0]
                ),
                "n_plus_components_available": int(
                    components[
                        components["component"].isin(PLUS_COMPONENTS)
                        & components["available"].astype(bool)
                    ].shape[0]
                ),
            }
        ]
    )
    observed.to_csv(paths.out_dir / "mbpi_observed_summary.csv", index=False)
    return components, classifier_components


def update_observed_summary_with_ci(paths: MbpiPaths) -> None:
    observed_path = paths.out_dir / "mbpi_observed_summary.csv"
    summary_path = paths.out_dir / "mbpi_summary.csv"
    if not observed_path.exists() or not summary_path.exists():
        return
    observed = pd.read_csv(observed_path)
    summary = pd.read_csv(summary_path)
    mbpi_row = summary.loc[summary["metric"] == "MBPI"].iloc[0]
    plus_row = summary.loc[summary["metric"] == "MBPI_plus"].iloc[0]
    observed["MBPI_ci_low"] = mbpi_row["q2.5"]
    observed["MBPI_ci_high"] = mbpi_row["q97.5"]
    observed["classification"] = mbpi_row["classification"]
    observed["MBPI_plus_ci_low"] = plus_row["q2.5"]
    observed["MBPI_plus_ci_high"] = plus_row["q97.5"]
    observed["classification_plus"] = plus_row["classification"]
    observed.to_csv(observed_path, index=False)


def derive_effect_vectors(df: pd.DataFrame, feature_cols: list[str]) -> dict[str, object]:
    X = df[feature_cols].to_numpy(float)
    y = df["PD"].to_numpy(int)
    cohort = df["Cohort"].to_numpy(str)
    d_full = build_design(y, cohort, include_pd=True, include_cohort=True)
    coef, *_ = np.linalg.lstsq(d_full, X, rcond=None)
    common = coef[-1].copy()
    k = min(10, len(feature_cols))
    top_idx = np.argsort(np.abs(common))[::-1][:k]
    sparse = np.zeros_like(common)
    sparse[top_idx] = common[top_idx]
    sparse = sparse - sparse.mean()
    common_norm = np.linalg.norm(sparse)
    if common_norm <= 0:
        common_norm = 1.0
    common_unit = sparse / common_norm

    cohort_vectors = {}
    neutral = X.copy()
    for level in sorted(pd.unique(cohort)):
        mask = cohort == level
        pd_mask = mask & (y == 1)
        ctrl_mask = mask & (y == 0)
        delta = X[pd_mask].mean(axis=0) - X[ctrl_mask].mean(axis=0)
        neutral[pd_mask] = neutral[pd_mask] - delta
        cohort_sparse = np.zeros_like(delta)
        top_c = np.argsort(np.abs(delta))[::-1][:k]
        cohort_sparse[top_c] = delta[top_c]
        cohort_sparse = cohort_sparse - cohort_sparse.mean()
        norm_c = np.linalg.norm(cohort_sparse)
        if norm_c <= 0:
            cohort_unit = common_unit.copy()
        else:
            cohort_unit = cohort_sparse / norm_c
        if np.dot(cohort_unit, common_unit) < 0:
            cohort_unit = -cohort_unit
        cohort_vectors[level] = cohort_unit
    neutral = neutral - neutral.mean(axis=1, keepdims=True)
    return {
        "neutral_X": neutral,
        "common_unit": common_unit,
        "common_norm": common_norm,
        "cohort_units": cohort_vectors,
        "top_features": [feature_cols[i] for i in top_idx],
    }


def stratified_bootstrap_indices(df: pd.DataFrame, rng: np.random.Generator) -> np.ndarray:
    indices = []
    for (_, _), group in df.groupby(["Cohort", "PD"], sort=True):
        group_idx = group.index.to_numpy()
        indices.append(rng.choice(group_idx, size=len(group_idx), replace=True))
    return np.concatenate(indices)


def component_rows_from_scores(scores: dict[str, float], domain_score: float = 0.5) -> pd.DataFrame:
    rows = []
    for component in PRIMARY_COMPONENTS:
        rows.append(_component_row(component, scores.get(component, np.nan), pd.notna(scores.get(component, np.nan)), "", "simulation", "", True, True))
    rows.append(_component_row("S_domain", domain_score, True, "", "simulation", "No correction procedure evaluated in simulation.", False, True))
    return pd.DataFrame(rows)


def compute_matrix_mbpi(X: np.ndarray,
                        y: Iterable[int],
                        cohort: Iterable[str],
                        feature_cols: list[str],
                        seed: int,
                        recompute_prediction: bool = True,
                        recompute_heterogeneity: bool = True,
                        fixed_scores: dict[str, float] | None = None) -> dict[str, float]:
    fixed_scores = fixed_scores or {}
    r2 = marginal_r2_cddr(X, y, cohort)
    geom = disease_axis_geometry(X, y, cohort)
    scores = {
        "S_CDDR": score_cddr(r2["CDDR"]),
        "S_angle": score_angle(geom["Angle_Degrees"]),
    }
    raw = {
        "R2_cohort": r2["R2_cohort"],
        "R2_disease": r2["R2_disease"],
        "CDDR": r2["CDDR"],
        "angle_degrees": geom["Angle_Degrees"],
        "pc1_variance_pct": geom["PC1_Variance_Explained_Pct"],
    }
    if recompute_prediction:
        pred = compute_prediction_components(X, y, cohort, feature_cols, seed=seed)
        scores.update({k: pred[k] for k in ["S_BPG", "S_LOCO", "S_feature"]})
        raw.update(pred)
    else:
        for k in ["S_BPG", "S_LOCO", "S_feature"]:
            scores[k] = fixed_scores.get(k, np.nan)
    if recompute_heterogeneity:
        het = heterogeneity_from_matrix(X, y, cohort)
        scores["S_heterogeneity"] = het["score"]
        raw.update(
            {
                "interaction_q05": het["interaction_count"],
                "interaction_n": het["interaction_n"],
                "cochran_q_q05": het["cochran_q_count"],
                "cochran_q_n": het["cochran_q_n"],
            }
        )
    else:
        scores["S_heterogeneity"] = fixed_scores.get("S_heterogeneity", np.nan)
    comp = component_rows_from_scores(scores, domain_score=fixed_scores.get("S_domain", 0.5))
    raw["MBPI"] = aggregate_scores(comp, include_domain=False)
    raw["MBPI_plus"] = aggregate_scores(comp, include_domain=True)
    raw.update(scores)
    return raw


def optimal_accuracy_at_threshold(y_true: np.ndarray, scores: np.ndarray, threshold: float = 0.66) -> float:
    pred = (scores > threshold).astype(int)
    return float(np.mean(pred == y_true))


def auprc_from_scores(y_true: np.ndarray, scores: np.ndarray) -> float:
    if len(np.unique(y_true)) < 2:
        return np.nan
    return float(average_precision_score(y_true, scores))


def roc_auc_from_scores(y_true: np.ndarray, scores: np.ndarray) -> float:
    if len(np.unique(y_true)) < 2:
        return np.nan
    return float(roc_auc_score(y_true, scores))


def command_string(argv: list[str] | None = None) -> str:
    argv = list(sys.argv if argv is None else argv)
    if argv and argv[0].endswith(".py"):
        argv[0] = argv[0].replace("\\", "/")
    return "python " + " ".join(argv)
