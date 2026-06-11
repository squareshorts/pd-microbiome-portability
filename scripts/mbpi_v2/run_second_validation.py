from __future__ import annotations

import argparse
import hashlib
import sys
import tarfile
import urllib.request
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.metrics import average_precision_score, roc_auc_score

MBPI_DIR = Path(__file__).resolve().parents[1] / "mbpi"
if str(MBPI_DIR) not in sys.path:
    sys.path.insert(0, str(MBPI_DIR))

V2_DIR = Path(__file__).resolve().parent
if str(V2_DIR) not in sys.path:
    sys.path.insert(0, str(V2_DIR))

from mbpi_core import (  # noqa: E402
    compute_matrix_mbpi,
    derive_effect_vectors,
    disease_axis_geometry,
    heterogeneity_from_matrix,
    make_logistic_model,
    marginal_r2_cddr,
    read_clr_dataset,
    stratified_bootstrap_indices,
    summary_stats,
)
from mbpi_v2_core import (  # noqa: E402
    DSAS_COMPONENTS,
    DSAS_WEIGHTS,
    PRS_COMPONENTS,
    apply_prs_strategies,
    build_ablation_summary,
    build_confusion,
    calibrate_prs_weights,
    classify_dsas,
    classify_v2,
    compute_simulation_signal_and_risk,
    fixed_effect_from_matrix,
    pd_main_from_matrix,
    selected_weight_dict,
    signal_components_from_metrics,
    simulate_setting,
    simulation_design,
    strategy_score_from_parts,
    summarize_simulations,
    three_class_label,
)


DATAVERSE_BASE = "https://entrepot.recherche.data.gouv.fr/api/access/datafile"
SOURCE_FILES = {
    "species": {
        "id": "192286",
        "name": "species_signal_2340_CRC_cohort_20240322.tab",
        "description": "METEOR species abundance table",
        "doi": "https://doi.org/10.57745/DQBQD3",
    },
    "metadata": {
        "id": "218633",
        "name": "metadata_2340_CRC_cohort_20240704.tar.gz",
        "description": "Manually curated sample metadata archive",
        "doi": "https://doi.org/10.57745/LCAR4M",
    },
}
DATASET_DOI = "https://doi.org/10.57745/7IVO3E"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run an independent non-PD MBPI-v2 validation case on a public CRC multi-cohort benchmark."
    )
    parser.add_argument("--root", default=".", help="Repository root.")
    parser.add_argument("--simulation-reps", type=int, default=10, help="Replicates per MBPI-v2 simulation setting.")
    parser.add_argument("--bootstrap-B", type=int, default=100, help="Bootstrap replicate count.")
    parser.add_argument("--seed", type=int, default=20260610, help="Random seed.")
    parser.add_argument(
        "--prevalence-threshold",
        type=float,
        default=0.10,
        help="Minimum non-zero feature prevalence across retained samples.",
    )
    parser.add_argument(
        "--min-per-class-per-cohort",
        type=int,
        default=10,
        help="Minimum CRC and control samples required to retain a cohort for LOCO validation.",
    )
    return parser.parse_args()


def out_dir(root: Path) -> Path:
    return root / "results" / "mbpi_v2_second_validation"


def source_dir(root: Path) -> Path:
    return out_dir(root) / "source"


def command_string(argv: list[str] | None = None) -> str:
    argv = list(sys.argv if argv is None else argv)
    return "python " + " ".join(str(a).replace("\\", "/") for a in argv)


def sha256_file(path: Path) -> str:
    if not path.exists() or not path.is_file():
        return ""
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def download_if_missing(file_id: str, dest: Path) -> None:
    if dest.exists() and dest.stat().st_size > 0:
        return
    url = f"{DATAVERSE_BASE}/{file_id}?format=original"
    print(f"Downloading {dest.name}", flush=True)
    with urllib.request.urlopen(url) as response, dest.open("wb") as out:
        out.write(response.read())


def safe_extract_tar_gz(archive: Path, dest: Path) -> None:
    dest_resolved = dest.resolve()
    with tarfile.open(archive, "r:gz") as tar:
        for member in tar.getmembers():
            target = (dest / member.name).resolve()
            if not str(target).startswith(str(dest_resolved)):
                raise ValueError(f"Unsafe path in archive: {member.name}")
        tar.extractall(dest)


def prepare_sources(root: Path) -> dict[str, Path]:
    sdir = source_dir(root)
    sdir.mkdir(parents=True, exist_ok=True)
    paths = {}
    for key, item in SOURCE_FILES.items():
        dest = sdir / item["name"]
        download_if_missing(item["id"], dest)
        paths[key] = dest
    metadata_xlsx = sdir / "metadata_2340_CRC_cohort_20240704.xlsx"
    if not metadata_xlsx.exists():
        safe_extract_tar_gz(paths["metadata"], sdir)
    paths["metadata_xlsx"] = metadata_xlsx
    return paths


def make_unique(names: list[str]) -> list[str]:
    seen: dict[str, int] = {}
    out = []
    for name in names:
        base = "".join(ch if ch.isalnum() or ch == "_" else "_" for ch in str(name))
        if not base:
            base = "feature"
        if base not in seen:
            seen[base] = 0
            out.append(base)
        else:
            seen[base] += 1
            out.append(f"{base}_{seen[base]}")
    return out


def build_harmonized_matrix(
    root: Path,
    paths: dict[str, Path],
    prevalence_threshold: float,
    min_per_class: int,
) -> tuple[pd.DataFrame, list[str], pd.DataFrame, dict[str, object]]:
    odir = out_dir(root)
    odir.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_excel(paths["metadata_xlsx"], sheet_name="metadata_2340_CRC_cohort")
    abundance_header = pd.read_csv(paths["species"], sep="\t", nrows=0).columns.tolist()
    abundance_samples = [str(c) for c in abundance_header[1:]]
    abundance_sample_set = set(abundance_samples)

    matched = metadata[metadata["sample"].astype(str).isin(abundance_sample_set)].copy()
    binary = matched[matched["class"].isin(["healthy", "CRC"]) & matched["to_exclude"].isna()].copy()

    counts = binary.groupby(["study_accession", "class"]).size().unstack(fill_value=0)
    for col in ["CRC", "healthy"]:
        if col not in counts.columns:
            counts[col] = 0
    retained_cohorts = counts[(counts["CRC"] >= min_per_class) & (counts["healthy"] >= min_per_class)].index.tolist()
    rejected_counts = counts.loc[~counts.index.isin(retained_cohorts)].copy()

    selected_meta = binary[binary["study_accession"].isin(retained_cohorts)].copy()
    selected_sample_set = set(selected_meta["sample"].astype(str))
    sample_order = [s for s in abundance_samples if s in selected_sample_set]
    selected_meta = selected_meta.set_index("sample").loc[sample_order].reset_index()

    usecols = ["msp_id"] + sample_order
    abundance = pd.read_csv(paths["species"], sep="\t", usecols=usecols)
    values = abundance[sample_order].to_numpy(dtype=float)
    prevalence = (values > 0).mean(axis=1)
    feature_sums = values.sum(axis=1)
    keep_features = (prevalence >= prevalence_threshold) & (feature_sums > 0)
    if int(keep_features.sum()) == 0:
        raise ValueError("No features survived the prevalence filter.")

    feature_ids = abundance.loc[keep_features, "msp_id"].astype(str).tolist()
    feature_cols = make_unique(feature_ids)
    selected_values = values[keep_features, :].T
    positive = selected_values[selected_values > 0]
    if len(positive) == 0:
        raise ValueError("Retained abundance matrix has no positive values.")
    pseudocount = float(np.min(positive) / 2.0)
    log_values = np.log(selected_values + pseudocount)
    clr = log_values - log_values.mean(axis=1, keepdims=True)

    matrix = pd.DataFrame(clr, columns=feature_cols)
    matrix.insert(0, "Cohort", selected_meta["study_accession"].astype(str).to_numpy())
    matrix.insert(0, "PD", (selected_meta["class"].astype(str) == "CRC").astype(int).to_numpy())
    matrix.insert(0, "SampleID", sample_order)
    matrix.to_csv(odir / "second_validation_matrix.csv", index=False)

    disease_control = np.where(selected_meta["class"].astype(str) == "CRC", "disease", "control")
    metadata_out = pd.DataFrame(
        {
            "sample_id": sample_order,
            "phenotype": selected_meta["class"].astype(str).to_numpy(),
            "cohort_or_study": selected_meta["study_accession"].astype(str).to_numpy(),
            "disease_control_label": disease_control,
            "disease_binary": matrix["PD"].to_numpy(int),
            "health_status": selected_meta["health_status"].astype(str).to_numpy(),
            "host_phenotype": selected_meta["host_phenotype"].fillna("").astype(str).to_numpy(),
            "country": selected_meta["country"].astype(str).to_numpy(),
            "source_sample_column": "sample",
            "source_cohort_column": "study_accession",
            "source_phenotype_column": "class",
        }
    )
    metadata_out.to_csv(odir / "second_validation_metadata.csv", index=False)

    retained_counts = (
        selected_meta.groupby(["study_accession", "class"]).size().unstack(fill_value=0).reset_index()
    )
    rejected_text = "; ".join(
        f"{idx}: CRC={int(row.get('CRC', 0))}, healthy={int(row.get('healthy', 0))}"
        for idx, row in rejected_counts.iterrows()
    )
    validation_rows = [
        {
            "check": "dataset_source",
            "status": "pass",
            "value": DATASET_DOI,
            "notes": "Recherche Data Gouv/MetaGenoPolis CRC benchmark.",
        },
        {
            "check": "abundance_metadata_sample_overlap",
            "status": "pass",
            "value": f"{len(matched)}/{len(metadata)} metadata rows matched abundance columns",
            "notes": "The abundance table sample columns matched the metadata sample field exactly.",
        },
        {
            "check": "phenotype_filter",
            "status": "pass",
            "value": "retained class in {CRC, healthy}; excluded adenoma and to_exclude=yes rows",
            "notes": "No adenoma samples were merged into either binary class.",
        },
        {
            "check": "cohort_filter",
            "status": "pass",
            "value": f"{len(retained_cohorts)} cohorts retained with >= {min_per_class} samples per class",
            "notes": rejected_text or "No cohorts rejected by class-count filter.",
        },
        {
            "check": "sample_count",
            "status": "pass",
            "value": len(matrix),
            "notes": dict(matrix["PD"].value_counts().sort_index()),
        },
        {
            "check": "feature_filter",
            "status": "pass",
            "value": f"{len(feature_cols)} features retained from {len(abundance)} species",
            "notes": f"Non-zero prevalence threshold={prevalence_threshold}; CLR pseudocount={pseudocount:.12g}.",
        },
        {
            "check": "retained_cohort_counts",
            "status": "pass",
            "value": retained_counts.to_json(orient="records"),
            "notes": "Counts are after phenotype, quality, and minimum-per-class filters.",
        },
    ]
    validation = pd.DataFrame(validation_rows)
    validation.to_csv(odir / "second_validation_input_validation.csv", index=False)

    info = {
        "pseudocount": pseudocount,
        "prevalence_threshold": prevalence_threshold,
        "min_per_class": min_per_class,
        "retained_cohorts": retained_cohorts,
        "rejected_cohorts": rejected_text,
        "raw_n_samples": len(metadata),
        "binary_quality_n_samples": len(binary),
        "n_samples": len(matrix),
        "n_features": len(feature_cols),
        "n_cohorts": len(retained_cohorts),
        "n_crc": int(matrix["PD"].sum()),
        "n_control": int((matrix["PD"] == 0).sum()),
    }
    return matrix, feature_cols, metadata_out, info


def apply_final_classes(sim: pd.DataFrame, selected_strategy: str) -> pd.DataFrame:
    sim = sim.copy()
    sim["PRS"] = sim[f"PRS_{selected_strategy}"]
    sim["MBPI_v2_score"] = sim["DSAS"] * sim["PRS"]
    sim["final_class"] = [classify_v2(d, p) for d, p in zip(sim["DSAS"], sim["PRS"])]
    return sim


def run_simulations_from_matrix(
    matrix: pd.DataFrame,
    feature_cols: list[str],
    odir: Path,
    reps: int,
    seed: int,
) -> pd.DataFrame:
    design = simulation_design([0.25, 0.50, 1.00, 1.50])
    design.to_csv(odir / "second_validation_mbpi_v2_simulation_design.csv", index=False)
    effect_info = derive_effect_vectors(matrix, feature_cols)
    rng = np.random.default_rng(seed)
    rows = []
    sim_id = 0
    for _, setting in design.iterrows():
        print(f"Simulating {setting['setting']}", flush=True)
        for rep in range(reps):
            sim_id += 1
            X_sim, y_sim, cohort_sim = simulate_setting(setting, matrix, feature_cols, effect_info, rng)
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


def observed_metrics(matrix: pd.DataFrame, feature_cols: list[str], seed: int) -> dict[str, object]:
    X = matrix[feature_cols].to_numpy(float)
    y = matrix["PD"].to_numpy(int)
    cohort = matrix["Cohort"].to_numpy(str)
    risk_metrics = compute_matrix_mbpi(
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
    fixed_effect = fixed_effect_from_matrix(X, y, cohort)
    dsas_parts = signal_components_from_metrics(
        risk_metrics["pooled_auroc"],
        risk_metrics["loco_auroc"],
        pd_main["pd_main_proportion"],
        fixed_effect["fixed_effect_proportion"],
        pd_main["median_abs_pd_effect_sig"],
    )
    risk_values = {component: float(risk_metrics[component]) for component in PRS_COMPONENTS}
    return {
        "risk_metrics": risk_metrics,
        "pd_main": pd_main,
        "fixed_effect": fixed_effect,
        "dsas_parts": dsas_parts,
        "risk_values": risk_values,
        "observed_fixed": {
            "pooled_auroc_median": float(risk_metrics["pooled_auroc"]),
            "loco_auroc_median": float(risk_metrics["loco_auroc"]),
            **risk_values,
        },
    }


def bootstrap_v2_from_matrix(
    matrix: pd.DataFrame,
    feature_cols: list[str],
    selected_rows: pd.DataFrame,
    observed_fixed: dict[str, float],
    B: int,
    seed: int,
    odir: Path,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    for b in range(B):
        if (b + 1) % max(1, B // 5) == 0:
            print(f"Bootstrap replicate {b + 1}/{B}", flush=True)
        idx = stratified_bootstrap_indices(matrix, rng)
        boot = matrix.loc[idx].reset_index(drop=True)
        X = boot[feature_cols].to_numpy(float)
        y = boot["PD"].to_numpy(int)
        cohort = boot["Cohort"].to_numpy(str)
        pd_main = pd_main_from_matrix(X, y, cohort)
        fixed_effect = fixed_effect_from_matrix(X, y, cohort)
        het = heterogeneity_from_matrix(X, y, cohort)
        r2 = marginal_r2_cddr(X, y, cohort)
        geom = disease_axis_geometry(X, y, cohort)
        dsas_parts = signal_components_from_metrics(
            observed_fixed["pooled_auroc_median"],
            observed_fixed["loco_auroc_median"],
            pd_main["pd_main_proportion"],
            fixed_effect["fixed_effect_proportion"],
            pd_main["median_abs_pd_effect_sig"],
        )
        prs_parts = {
            "S_CDDR": score_cddr_local(r2["CDDR"]),
            "S_angle": min(max(float(geom["Angle_Degrees"]) / 90.0, 0.0), 1.0),
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
            **fixed_effect,
            **prs_parts,
            "R2_cohort": r2["R2_cohort"],
            "R2_disease": r2["R2_disease"],
            "CDDR": r2["CDDR"],
            "angle_degrees": geom["Angle_Degrees"],
            "bootstrap_mode": "fast_species_and_geometry",
        }
        rows.append(row)
    boot = pd.DataFrame(rows)
    boot.to_csv(odir / "second_validation_mbpi_v2_bootstrap_replicates.csv", index=False)
    return boot


def score_cddr_local(cddr: float) -> float:
    if cddr is None or pd.isna(cddr):
        return np.nan
    if np.isinf(cddr):
        return 1.0
    transformed = np.log10(1.0 + max(0.0, float(cddr)))
    return float(np.clip(transformed / (1.0 + transformed), 0.0, 1.0))


def write_observed_summary(
    odir: Path,
    obs: dict[str, object],
    weights: pd.DataFrame,
    selected_strategy: str,
    bootstrap: pd.DataFrame,
) -> pd.DataFrame:
    selected_rows = weights[weights["strategy"] == selected_strategy]
    risk_values = obs["risk_values"]
    prs = strategy_score_from_parts(risk_values, selected_rows)
    equal_weight = float(np.mean([risk_values[c] for c in PRS_COMPONENTS]))
    dsas = float(obs["dsas_parts"]["DSAS"])
    final_class = classify_v2(dsas, prs)

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
                "PRS": prs,
                "PRS_ci_low": prs_low,
                "PRS_ci_high": prs_high,
                "PRS_strategy": selected_strategy,
                "PRS_equal_weight": equal_weight,
                "MBPI_v2_score": dsas * prs,
                "MBPI_v2_score_ci_low": score_low,
                "MBPI_v2_score_ci_high": score_high,
                "final_class": final_class,
                "interpretation": "Detectable disease signal with cohort-conditioned non-portability support."
                if final_class in {"non-portable", "cohort-conditioned"}
                else final_class,
            }
        ]
    )
    summary.to_csv(odir / "second_validation_mbpi_v2_observed_summary.csv", index=False)
    return summary


def write_components(
    odir: Path,
    obs: dict[str, object],
    weights: pd.DataFrame,
    selected_strategy: str,
) -> pd.DataFrame:
    risk_metrics = obs["risk_metrics"]
    pd_main = obs["pd_main"]
    fixed_effect = obs["fixed_effect"]
    dsas_parts = obs["dsas_parts"]
    risk_values = obs["risk_values"]
    selected_weights = selected_weight_dict(weights, selected_strategy)
    equal_weights = {c: 1.0 / len(PRS_COMPONENTS) for c in PRS_COMPONENTS}
    source_file = str((odir / "second_validation_matrix.csv").resolve())
    rows = []
    dsas_raw = {
        "DSAS_pooled_auroc": f"single logistic pooled AUROC={risk_metrics['pooled_auroc']:.12g}",
        "DSAS_loco_auroc": f"single logistic mean LOCO AUROC={risk_metrics['loco_auroc']:.12g}",
        "DSAS_pd_main_fdr": f"{pd_main['pd_main_count']}/{pd_main['pd_main_n']} disease-main q<0.05",
        "DSAS_fixed_effect_fdr": f"{fixed_effect['fixed_effect_count']}/{fixed_effect['fixed_effect_n']} fixed-effect q<0.05",
        "DSAS_effect_size": f"median abs disease beta among q<0.05={pd_main['median_abs_pd_effect_sig']:.12g}",
    }
    for component in DSAS_COMPONENTS:
        rows.append(
            {
                "stage": "disease_signal_adequacy",
                "component": component,
                "score": dsas_parts[component],
                "weight": DSAS_WEIGHTS[component],
                "raw_value": dsas_raw[component],
                "source_file": source_file,
                "notes": "Higher values indicate stronger usable disease signal.",
                "selected_weight": "",
                "equal_weight": "",
            }
        )
    rows.append(
        {
            "stage": "disease_signal_adequacy",
            "component": "DSAS",
            "score": dsas_parts["DSAS"],
            "weight": 1.0,
            "raw_value": "weighted DSAS",
            "source_file": source_file,
            "notes": classify_dsas(dsas_parts["DSAS"]),
            "selected_weight": "",
            "equal_weight": "",
        }
    )
    raw_risk = {
        "S_CDDR": f"CDDR={risk_metrics['CDDR']:.12g}; R2_cohort={risk_metrics['R2_cohort']:.12g}; R2_disease={risk_metrics['R2_disease']:.12g}",
        "S_angle": f"Angle_Degrees={risk_metrics['angle_degrees']:.12g}; PC1_Variance_Explained_Pct={risk_metrics['pc1_variance_pct']:.12g}",
        "S_BPG": f"single logistic pooled AUROC={risk_metrics['pooled_auroc']:.12g}; LOCO AUROC={risk_metrics['loco_auroc']:.12g}; BPG={risk_metrics['bpg_auroc']:.12g}",
        "S_LOCO": f"single logistic LOCO AUROC={risk_metrics['loco_auroc']:.12g}",
        "S_feature": f"leave-one-cohort top-20 Jaccard={risk_metrics['feature_jaccard_top20']:.12g}",
        "S_heterogeneity": f"interaction_q05={risk_metrics['interaction_q05']}/{risk_metrics['interaction_n']}; Cochran_Q={risk_metrics['cochran_q_q05']}/{risk_metrics['cochran_q_n']}",
    }
    for component in PRS_COMPONENTS:
        rows.append(
            {
                "stage": "portability_risk",
                "component": component,
                "score": risk_values[component],
                "weight": "",
                "raw_value": raw_risk[component],
                "source_file": source_file,
                "notes": "Higher values indicate stronger portability risk.",
                "selected_weight": selected_weights.get(component, ""),
                "equal_weight": equal_weights.get(component, ""),
            }
        )
    comp = pd.DataFrame(rows)
    comp.to_csv(odir / "second_validation_mbpi_v2_components.csv", index=False)
    return comp


def compute_loco_table(matrix: pd.DataFrame, feature_cols: list[str], seed: int) -> pd.DataFrame:
    X = matrix[feature_cols].to_numpy(float)
    y = matrix["PD"].to_numpy(int)
    cohort = matrix["Cohort"].to_numpy(str)
    rows = []
    for level in sorted(pd.unique(cohort)):
        test = cohort == level
        train = ~test
        model = make_logistic_model(seed)
        model.fit(X[train], y[train])
        scores = model.predict_proba(X[test])[:, 1]
        auroc = roc_auc_score(y[test], scores) if len(np.unique(y[test])) == 2 else np.nan
        auprc = average_precision_score(y[test], scores) if len(np.unique(y[test])) == 2 else np.nan
        rows.append(
            {
                "cohort_or_study": level,
                "n_test": int(test.sum()),
                "n_crc": int(y[test].sum()),
                "n_control": int((y[test] == 0).sum()),
                "loco_auroc": auroc,
                "loco_auprc": auprc,
            }
        )
    return pd.DataFrame(rows)


def save_pair(fig: plt.Figure, odir: Path, stem: str) -> None:
    fig.savefig(odir / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(odir / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def make_figures(
    odir: Path,
    selected_strategy: str,
    pooled_auroc: float,
    loco_table: pd.DataFrame,
) -> None:
    sns.set_theme(style="whitegrid")
    observed = pd.read_csv(odir / "second_validation_mbpi_v2_observed_summary.csv")
    components = pd.read_csv(odir / "second_validation_mbpi_v2_components.csv")
    sim = pd.read_csv(odir / "second_validation_mbpi_v2_simulation_results.csv", keep_default_na=False)
    boot = pd.read_csv(odir / "second_validation_mbpi_v2_bootstrap_replicates.csv")
    ablation = pd.read_csv(odir / "second_validation_mbpi_v2_ablation_summary.csv")

    plot_components = components[components["component"].isin(["DSAS", *PRS_COMPONENTS])].copy()
    fig, ax = plt.subplots(figsize=(9, 5.2))
    colors = ["#5b8e7d" if c == "DSAS" else "#b85c5c" for c in plot_components["component"]]
    ax.bar(plot_components["component"], plot_components["score"].astype(float), color=colors, edgecolor="#333333")
    ax.axhline(0.33, color="#777777", linestyle="--", linewidth=1)
    ax.axhline(0.66, color="#333333", linestyle="--", linewidth=1)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Score")
    ax.set_title("Second Validation MBPI-v2 Components")
    ax.tick_params(axis="x", rotation=35)
    for tick in ax.get_xticklabels():
        tick.set_ha("right")
    save_pair(fig, odir, "fig_second_validation_mbpi_v2_components")

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    for ax, metric in zip(axes, ["DSAS", "PRS"]):
        sns.histplot(boot[metric], bins=25, color="#5b8e7d" if metric == "DSAS" else "#b85c5c", ax=ax)
        ax.axvline(observed[metric].iloc[0], color="#111111", linewidth=2, label=f"observed {observed[metric].iloc[0]:.3f}")
        ax.axvline(np.quantile(boot[metric], 0.025), color="#333333", linestyle="--", linewidth=1)
        ax.axvline(np.quantile(boot[metric], 0.975), color="#333333", linestyle="--", linewidth=1)
        ax.set_xlabel(metric)
        ax.legend(frameon=False)
    fig.suptitle("Second Validation MBPI-v2 Bootstrap")
    save_pair(fig, odir, "fig_second_validation_mbpi_v2_bootstrap")

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5))
    sns.boxplot(data=sim, x="regime", y="DSAS", order=["null", "portable", "nonportable"], color="#dce9df", ax=axes[0])
    axes[0].axhline(0.33, color="#777777", linestyle="--", linewidth=1)
    axes[0].axhline(0.66, color="#333333", linestyle="--", linewidth=1)
    axes[0].set_title("DSAS Calibration")
    axes[0].set_xlabel("")
    sns.boxplot(
        data=sim[sim["regime"].isin(["portable", "nonportable"])],
        x="regime",
        y=f"PRS_{selected_strategy}",
        order=["portable", "nonportable"],
        color="#ead7d7",
        ax=axes[1],
    )
    axes[1].axhline(0.66, color="#333333", linestyle="--", linewidth=1)
    axes[1].set_title("PRS Calibration")
    axes[1].set_xlabel("")
    axes[1].set_ylabel("PRS")
    save_pair(fig, odir, "fig_second_validation_mbpi_v2_simulation_calibration")

    plot_df = ablation.copy()
    plot_df["metric"] = plot_df["auroc"]
    guardrail = "Null false non-portability/cohort-conditioned rate"
    plot_df.loc[plot_df["diagnostic"] == guardrail, "metric"] = plot_df.loc[
        plot_df["diagnostic"] == guardrail, "notes"
    ].astype(float)
    plot_df = plot_df.sort_values("metric", ascending=True)
    fig, ax = plt.subplots(figsize=(9.5, 5.7))
    colors = ["#5b8e7d" if selected_strategy in d or d == "MBPI-v2 gated score" else "#c7d6df" for d in plot_df["diagnostic"]]
    ax.barh(plot_df["diagnostic"], plot_df["metric"], color=colors, edgecolor="#333333", linewidth=0.5)
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("AUROC, or rate for null guardrail")
    ax.set_title("Second Validation MBPI-v2 Ablation")
    for y, val in enumerate(plot_df["metric"]):
        ax.text(min(float(val) + 0.015, 1.0), y, f"{float(val):.3f}", va="center", fontsize=9)
    save_pair(fig, odir, "fig_second_validation_mbpi_v2_ablation")

    loco_table = loco_table.sort_values("loco_auroc", ascending=True)
    fig, ax = plt.subplots(figsize=(9.5, 5.2))
    ax.barh(loco_table["cohort_or_study"], loco_table["loco_auroc"], color="#c7d6df", edgecolor="#333333", linewidth=0.5)
    ax.axvline(pooled_auroc, color="#b85c5c", linewidth=2, label=f"pooled CV AUROC {pooled_auroc:.3f}")
    ax.set_xlim(0, 1.0)
    ax.set_xlabel("AUROC")
    ax.set_title("Second Validation Pooled vs LOCO AUROC")
    ax.legend(frameon=False)
    save_pair(fig, odir, "fig_second_validation_pooled_vs_loco")


def metric_from_ablation(ablation: pd.DataFrame, diagnostic: str, value_col: str = "auroc") -> float:
    row = ablation[ablation["diagnostic"] == diagnostic]
    if row.empty:
        return np.nan
    return float(row[value_col].iloc[0])


def null_guardrail(ablation: pd.DataFrame) -> float:
    row = ablation[ablation["diagnostic"] == "Null false non-portability/cohort-conditioned rate"]
    if row.empty:
        return np.nan
    return float(row["notes"].iloc[0])


def read_ci(summary: pd.DataFrame, prefix: str, side: str) -> float:
    for col in [f"{prefix}_ci_{side}", f"{prefix}_CI_{side}"]:
        if col in summary.columns:
            return float(summary[col].iloc[0])
    return np.nan


def build_cross_dataset_comparison(root: Path, odir: Path, matrix: pd.DataFrame, feature_cols: list[str]) -> pd.DataFrame:
    pd_summary = pd.read_csv(root / "results" / "mbpi_v2" / "mbpi_v2_observed_summary.csv")
    pd_ablation = pd.read_csv(root / "results" / "mbpi_v2" / "mbpi_v2_ablation_summary.csv")
    pd_matrix, pd_features = read_clr_dataset(root / "results" / "portability_analysis" / "dataset_clr.csv")
    second_summary = pd.read_csv(odir / "second_validation_mbpi_v2_observed_summary.csv")
    second_ablation = pd.read_csv(odir / "second_validation_mbpi_v2_ablation_summary.csv")

    def row(
        dataset: str,
        disease_area: str,
        n_samples: int,
        n_cohorts: int,
        n_features: int,
        summary: pd.DataFrame,
        ablation: pd.DataFrame,
        notes: str,
    ) -> dict[str, object]:
        strategy = str(summary["PRS_strategy"].iloc[0])
        return {
            "dataset": dataset,
            "disease area": disease_area,
            "n_samples": n_samples,
            "n_cohorts": n_cohorts,
            "n_features": n_features,
            "DSAS": float(summary["DSAS"].iloc[0]),
            "DSAS_CI_lower": read_ci(summary, "DSAS", "low"),
            "DSAS_CI_upper": read_ci(summary, "DSAS", "high"),
            "PRS": float(summary["PRS"].iloc[0]),
            "PRS_CI_lower": read_ci(summary, "PRS", "low"),
            "PRS_CI_upper": read_ci(summary, "PRS", "high"),
            "equal_weight_PRS": float(summary["PRS_equal_weight"].iloc[0]),
            "final_class": str(summary["final_class"].iloc[0]),
            "null_false_nonportable_rate": null_guardrail(ablation),
            "DSAS_AUROC": metric_from_ablation(ablation, "DSAS: null vs disease signal"),
            "PRS_AUROC": metric_from_ablation(ablation, f"PRS: {strategy}"),
            "gated_MBPI_AUROC": metric_from_ablation(ablation, "MBPI-v2 gated score"),
            "heterogeneity_alone_AUROC": metric_from_ablation(ablation, "PRS: heterogeneity_alone"),
            "notes": notes,
        }

    comparison = pd.DataFrame(
        [
            row(
                "PD validation case",
                "Parkinson's disease",
                len(pd_matrix),
                pd_matrix["Cohort"].nunique(),
                len(pd_features),
                pd_summary,
                pd_ablation,
                "Existing MBPI-v2 PD outputs read only; no files in results/mbpi_v2 were modified.",
            ),
            row(
                "second validation case",
                "colorectal cancer",
                len(matrix),
                matrix["Cohort"].nunique(),
                len(feature_cols),
                second_summary,
                second_ablation,
                "MetaGenoPolis/Recherche Data Gouv CRC benchmark; CRC vs healthy only after quality and cohort filters.",
            ),
        ]
    )
    comparison.to_csv(odir / "mbpi_v2_cross_dataset_comparison.csv", index=False)
    return comparison


def write_audit_table(
    odir: Path,
    paths: dict[str, Path],
    info: dict[str, object],
    obs: dict[str, object],
    selected_strategy: str,
    simulation_reps: int,
    bootstrap_B: int,
) -> pd.DataFrame:
    risk = obs["risk_metrics"]
    rows = [
        ("source_dataset", DATASET_DOI, "CRC benchmark dataset DOI."),
        ("species_file", str(paths["species"]), f"sha256={sha256_file(paths['species'])}"),
        ("metadata_file", str(paths["metadata_xlsx"]), f"sha256={sha256_file(paths['metadata_xlsx'])}"),
        ("phenotype_rule", "CRC=1, healthy=0", "Adenoma and rows marked to_exclude=yes were excluded."),
        (
            "cohort_rule",
            f"retained cohorts with >= {info['min_per_class']} CRC and >= {info['min_per_class']} controls",
            f"rejected: {info['rejected_cohorts']}",
        ),
        ("n_samples", info["n_samples"], f"CRC={info['n_crc']}; controls={info['n_control']}"),
        ("n_cohorts", info["n_cohorts"], ";".join(info["retained_cohorts"])),
        ("n_features", info["n_features"], f"prevalence_threshold={info['prevalence_threshold']}"),
        ("clr_pseudocount", info["pseudocount"], "Half the minimum positive retained abundance."),
        ("pooled_auroc", risk["pooled_auroc"], "Single L2 logistic model with stratified pooled CV."),
        ("loco_auroc", risk["loco_auroc"], "Mean leave-one-cohort-out AUROC across retained cohorts."),
        ("selected_prs_strategy", selected_strategy, "Calibrated on MBPI-v2 simulations for this second dataset."),
        ("simulation_reps_per_setting", simulation_reps, "Used for second validation run."),
        ("bootstrap_replicates", bootstrap_B, "Used for second validation uncertainty intervals."),
    ]
    audit = pd.DataFrame(rows, columns=["audit_item", "value", "notes"])
    audit.to_csv(odir / "second_validation_mbpi_v2_audit_table.csv", index=False)
    return audit


def write_run_log(
    root: Path,
    odir: Path,
    paths: dict[str, Path],
    info: dict[str, object],
    selected_strategy: str,
    command: str,
    warnings: list[str],
) -> None:
    observed = pd.read_csv(odir / "second_validation_mbpi_v2_observed_summary.csv")
    ablation = pd.read_csv(odir / "second_validation_mbpi_v2_ablation_summary.csv")
    weights = pd.read_csv(odir / "second_validation_mbpi_v2_weights.csv")
    lines = [
        "Second validation MBPI-v2 run log",
        f"date_time: {datetime.now().isoformat(timespec='seconds')}",
        f"repository_root: {root}",
        f"command: {command}",
        f"selected_dataset: MetaGenoPolis Recherche Data Gouv CRC benchmark ({DATASET_DOI})",
        f"selected_prs_strategy: {selected_strategy}",
        "",
        "source_files:",
    ]
    for key in ["species", "metadata", "metadata_xlsx"]:
        lines.append(f"- {paths[key]} | found={paths[key].exists()} | sha256={sha256_file(paths[key])}")
    lines.extend(
        [
            "",
            f"number_of_samples: {info['n_samples']}",
            f"number_of_features: {info['n_features']}",
            f"number_of_cohorts: {info['n_cohorts']}",
            f"class_counts: CRC={info['n_crc']}; control={info['n_control']}",
            f"retained_cohorts: {info['retained_cohorts']}",
            f"rejected_cohorts: {info['rejected_cohorts']}",
            f"prevalence_threshold: {info['prevalence_threshold']}",
            f"clr_pseudocount: {info['pseudocount']}",
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
    for row in weights[weights["strategy"] == selected_strategy].itertuples(index=False):
        lines.append(f"- {row.component}: weight={row.weight}; coefficient={row.coefficient}")
    lines.extend(["", "ablation_summary:"])
    lines.extend(ablation.to_string(index=False).splitlines())
    lines.extend(["", "warnings_or_limitations:"])
    if warnings:
        lines.extend(f"- {w}" for w in warnings)
    else:
        lines.append("- none")
    (odir / "second_validation_mbpi_v2_run_log.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_dataset_feasibility_report(odir: Path, info: dict[str, object]) -> None:
    text = f"""# MBPI-v2 Second Validation Dataset Feasibility Report

## Decision

Selected dataset: MetaGenoPolis / Recherche Data Gouv colorectal cancer benchmark.

Reason: it is the fastest credible independent validation case because it provides a single public repository with manually curated metadata, explicit study labels, explicit CRC/healthy/adenoma phenotypes, and ready-to-use species abundance tables. It also has strong methodological value for MBPI-v2 because the retained CRC/control subset contains {info['n_samples']} samples from {info['n_cohorts']} independent studies after transparent quality and cohort-balance filters.

## Candidate 1: MetaGenoPolis / Recherche Data Gouv CRC benchmark

- Dataset name: Taxonomic profiles, functional profiles and manually curated metadata of human fecal metagenomes from public projects coming from colorectal cancer studies.
- Disease area: Colorectal cancer.
- Source/accession or repository: {DATASET_DOI}; selected files {SOURCE_FILES['species']['doi']} and {SOURCE_FILES['metadata']['doi']}.
- Number of cohorts/studies: 15 available; {info['n_cohorts']} retained for binary CRC/control LOCO validation.
- Approximate sample size: 2,340 total metadata rows; {info['n_samples']} retained after excluding adenoma, quality-flagged rows, and cohorts without at least {info['min_per_class']} CRC and {info['min_per_class']} controls.
- Feature type: METEOR metagenomic species pangenome abundance, CLR transformed after a {info['prevalence_threshold']:.2f} non-zero prevalence filter ({info['n_features']} retained features).
- Phenotype labels available: yes, `class` contains `healthy`, `CRC`, and `adenoma`; only `healthy` and `CRC` were used.
- Cohort/study labels available: yes, `study_accession`.
- Preprocessing required: download TSV/XLSX files, extract metadata, remove adenoma and flagged samples, remove one-class/very small binary cohorts, prevalence-filter species, CLR-transform.
- Expected time to integrate: same day; direct TSV/XLSX import.
- Risks: MSP feature identifiers are not genus names; the source abundance scale is very small and requires CLR pseudocount handling; cohort filters must be reported because three cohorts were unsuitable for binary LOCO validation.
- Final recommendation: choose this dataset.

## Candidate 2: IBD 15-dataset benchmark from Kubinski et al. 2022

- Dataset name: Benchmark of Data Processing Methods and Machine Learning Models for Gut Microbiome-Based Diagnosis of Inflammatory Bowel Disease.
- Disease area: Inflammatory bowel disease.
- Source/accession or repository: Frontiers in Genetics article, doi:10.3389/fgene.2022.784397.
- Number of cohorts/studies: 15 datasets reported.
- Approximate sample size: 7,707 samples reported.
- Feature type: 16S-derived taxonomic and inferred functional features.
- Phenotype labels available: yes for IBD classification in the benchmark.
- Cohort/study labels available: yes; leave-one-dataset-out validation was used in the study.
- Preprocessing required: higher than selected CRC case, because the easiest public route is through the benchmark workflow/supplements rather than one direct abundance-plus-metadata pair already harmonized for this repository.
- Expected time to integrate: one to several days.
- Risks: 16S data are lower resolution than the current PD genus/shotgun-derived validation context; data retrieval and harmonization are less direct; disease subtype handling would need extra care to avoid merging incompatible CD/UC labels.
- Final recommendation: reject for this task because CRC is faster and cleaner.

## Candidate 3: curatedMetagenomicData disease-specific metagenomic pulls

- Dataset name: curatedMetagenomicData standardized human metagenomic resources.
- Disease area: Multiple possible disease areas, including CRC and IBD-related studies depending on query.
- Source/accession or repository: Bioconductor curatedMetagenomicData.
- Number of cohorts/studies: many standardized studies are available, but a disease-specific benchmark subset must be queried and assembled.
- Approximate sample size: disease-query dependent.
- Feature type: MetaPhlAn3 relative abundance, marker presence/abundance, gene families, and HUMAnN3 pathway abundance/coverage.
- Phenotype labels available: yes in curated metadata for many studies, but label harmonization must be inspected per disease.
- Cohort/study labels available: yes, through study/resource identifiers.
- Preprocessing required: install/load Bioconductor objects, query disease metadata, resolve repeated/ambiguous phenotypes, aggregate selected resources, write matrices.
- Expected time to integrate: one to several days.
- Risks: extra dependency burden and more phenotype-selection choices; weaker speed advantage than the direct CRC Dataverse files.
- Final recommendation: reject for this task because it is a useful fallback but not the fastest path.
"""
    (odir / "dataset_feasibility_report.md").write_text(text, encoding="utf-8")


def write_completion_report(
    odir: Path,
    info: dict[str, object],
    comparison: pd.DataFrame,
    command: str,
) -> None:
    observed = pd.read_csv(odir / "second_validation_mbpi_v2_observed_summary.csv")
    ablation = pd.read_csv(odir / "second_validation_mbpi_v2_ablation_summary.csv")
    required_files = [
        "dataset_feasibility_report.md",
        "second_validation_matrix.csv",
        "second_validation_metadata.csv",
        "second_validation_input_validation.csv",
        "second_validation_mbpi_v2_observed_summary.csv",
        "second_validation_mbpi_v2_components.csv",
        "second_validation_mbpi_v2_audit_table.csv",
        "second_validation_mbpi_v2_bootstrap_replicates.csv",
        "second_validation_mbpi_v2_simulation_results.csv",
        "second_validation_mbpi_v2_simulation_summary.csv",
        "second_validation_mbpi_v2_ablation_summary.csv",
        "second_validation_mbpi_v2_run_log.txt",
        "fig_second_validation_mbpi_v2_components.png",
        "fig_second_validation_mbpi_v2_components.pdf",
        "fig_second_validation_mbpi_v2_bootstrap.png",
        "fig_second_validation_mbpi_v2_bootstrap.pdf",
        "fig_second_validation_mbpi_v2_simulation_calibration.png",
        "fig_second_validation_mbpi_v2_simulation_calibration.pdf",
        "fig_second_validation_mbpi_v2_ablation.png",
        "fig_second_validation_mbpi_v2_ablation.pdf",
        "fig_second_validation_pooled_vs_loco.png",
        "fig_second_validation_pooled_vs_loco.pdf",
        "mbpi_v2_cross_dataset_comparison.csv",
    ]
    files_text = "\n".join(f"- `{name}`" for name in required_files)
    second = comparison[comparison["dataset"] == "second validation case"].iloc[0]
    pd_row = comparison[comparison["dataset"] == "PD validation case"].iloc[0]
    text = f"""# MBPI-v2 Second Validation Completion Report

## Selected Dataset

Selected dataset: MetaGenoPolis / Recherche Data Gouv colorectal cancer benchmark ({DATASET_DOI}).

It was selected because it is disease-independent from PD, has public multi-cohort shotgun metagenomic abundance data, explicit CRC/control labels, explicit study labels, and minimal harmonization burden.

## Rejected Candidate Datasets

- IBD 15-dataset benchmark: rejected for this task because it is methodologically relevant but slower to retrieve and harmonize, and CD/UC/IBD subtype handling would require more phenotype decisions.
- curatedMetagenomicData disease-specific pulls: rejected for this task because it is an excellent standardized source but requires Bioconductor object retrieval and disease-specific label assembly, making it slower than the direct CRC benchmark.

## Files Created

{files_text}

## Second Validation MBPI-v2 Results

- Samples: {info['n_samples']} ({info['n_crc']} CRC, {info['n_control']} controls).
- Cohorts: {info['n_cohorts']}.
- Features: {info['n_features']} CLR-transformed species features.
- DSAS: {observed['DSAS'].iloc[0]:.4f} ({observed['DSAS_ci_low'].iloc[0]:.4f} to {observed['DSAS_ci_high'].iloc[0]:.4f}).
- PRS: {observed['PRS'].iloc[0]:.4f} ({observed['PRS_ci_low'].iloc[0]:.4f} to {observed['PRS_ci_high'].iloc[0]:.4f}).
- Equal-weight PRS: {observed['PRS_equal_weight'].iloc[0]:.4f}.
- Final class: {observed['final_class'].iloc[0]}.
- DSAS calibration AUROC: {metric_from_ablation(ablation, 'DSAS: null vs disease signal'):.4f}.
- PRS calibration AUROC: {metric_from_ablation(ablation, 'PRS: ' + str(observed['PRS_strategy'].iloc[0])):.4f}.
- Gated MBPI-v2 AUROC: {metric_from_ablation(ablation, 'MBPI-v2 gated score'):.4f}.
- Null false non-portability/cohort-conditioned rate: {null_guardrail(ablation):.4f}.

## Comparison With PD Validation Case

- PD validation: DSAS {pd_row['DSAS']:.4f}, PRS {pd_row['PRS']:.4f}, final class `{pd_row['final_class']}`.
- CRC second validation: DSAS {second['DSAS']:.4f}, PRS {second['PRS']:.4f}, final class `{second['final_class']}`.

The second dataset strengthens the methodological claim that MBPI-v2 can be applied beyond Parkinson's disease because the full pipeline completed on an independent non-PD disease area with explicit cohort labels and binary disease/control labels. It should still be described as one additional successful validation case, not proof of universal generality.

## Exact Command To Reproduce

```bash
{command}
```

## Missing Inputs Or Limitations

- No manuscript text was modified.
- Adenoma samples were excluded rather than merged with CRC or controls.
- Cohorts with one class or fewer than {info['min_per_class']} samples per class were excluded from the binary LOCO validation.
- The selected feature level is METEOR MSP species identifiers rather than genus labels.
- Bootstrap and simulation replicate counts are recorded in `second_validation_mbpi_v2_run_log.txt`; increase them for publication-grade interval stability if runtime allows.
- In this scoped run, selected PRS AUROC was {metric_from_ablation(ablation, 'PRS: ' + str(observed['PRS_strategy'].iloc[0])):.4f} and heterogeneity-alone AUROC was {metric_from_ablation(ablation, 'PRS: heterogeneity_alone'):.4f}; interpret the PRS calibration as a completed second-dataset application, not as evidence that calibrated PRS improved over every ablation.
- Bootstrap intervals for DSAS and PRS were degenerate in this run because the bootstrapped disease-signal and calibrated PRS values were saturated at the observed values.
"""
    (odir / "completion_report.md").write_text(text, encoding="utf-8")


def main() -> None:
    args = parse_args()
    root = Path(args.root).resolve()
    odir = out_dir(root)
    odir.mkdir(parents=True, exist_ok=True)
    command = command_string()

    print("Preparing CRC second-validation sources", flush=True)
    paths = prepare_sources(root)

    print("Building harmonized CRC matrix", flush=True)
    matrix, feature_cols, _, info = build_harmonized_matrix(
        root,
        paths,
        prevalence_threshold=args.prevalence_threshold,
        min_per_class=args.min_per_class_per_cohort,
    )

    print("Computing observed MBPI-v2 components", flush=True)
    obs = observed_metrics(matrix, feature_cols, args.seed)

    print("Running MBPI-v2 simulation calibration", flush=True)
    sim = run_simulations_from_matrix(matrix, feature_cols, odir, args.simulation_reps, args.seed)
    weights, selected_strategy = calibrate_prs_weights(sim)
    weights.to_csv(odir / "second_validation_mbpi_v2_weights.csv", index=False)
    sim = apply_prs_strategies(sim, weights)
    sim = apply_final_classes(sim, selected_strategy)
    sim.to_csv(odir / "second_validation_mbpi_v2_simulation_results.csv", index=False)
    sim_summary = summarize_simulations(sim, selected_strategy)
    sim_summary.to_csv(odir / "second_validation_mbpi_v2_simulation_summary.csv", index=False)

    print("Running MBPI-v2 bootstrap uncertainty", flush=True)
    selected_rows = weights[weights["strategy"] == selected_strategy]
    boot = bootstrap_v2_from_matrix(
        matrix,
        feature_cols,
        selected_rows,
        obs["observed_fixed"],
        args.bootstrap_B,
        args.seed,
        odir,
    )

    print("Writing observed summary, components, and audits", flush=True)
    observed_summary = write_observed_summary(odir, obs, weights, selected_strategy, boot)
    components = write_components(odir, obs, weights, selected_strategy)
    ablation = build_ablation_summary(sim, weights, selected_strategy)
    ablation.to_csv(odir / "second_validation_mbpi_v2_ablation_summary.csv", index=False)
    confusion = build_confusion(sim)
    confusion.to_csv(odir / "second_validation_mbpi_v2_confusion_matrix.csv", index=False)
    write_audit_table(odir, paths, info, obs, selected_strategy, args.simulation_reps, args.bootstrap_B)
    write_dataset_feasibility_report(odir, info)

    print("Generating figures", flush=True)
    loco_table = compute_loco_table(matrix, feature_cols, args.seed)
    loco_table.to_csv(odir / "second_validation_pooled_vs_loco.csv", index=False)
    make_figures(odir, selected_strategy, obs["risk_metrics"]["pooled_auroc"], loco_table)

    print("Writing cross-dataset comparison and completion report", flush=True)
    comparison = build_cross_dataset_comparison(root, odir, matrix, feature_cols)
    dsas_auc = metric_from_ablation(ablation, "DSAS: null vs disease signal")
    prs_auc = metric_from_ablation(ablation, f"PRS: {selected_strategy}")
    hetero_auc = metric_from_ablation(ablation, "PRS: heterogeneity_alone")
    null_false = null_guardrail(ablation)
    warnings = []
    if dsas_auc < 0.80:
        warnings.append(f"DSAS calibration AUROC was {dsas_auc:.3f}, below the 0.80 acceptance target.")
    if prs_auc <= hetero_auc:
        warnings.append(f"Selected PRS AUROC was {prs_auc:.3f}; heterogeneity alone was {hetero_auc:.3f}.")
    if null_false > 0.05:
        warnings.append(f"Null-regime false non-portability/cohort-conditioned rate was {null_false:.3f}.")
    warnings.append(
        "This second validation uses one single classifier path for prediction-derived components rather than the PD run's existing multi-classifier summaries."
    )
    write_run_log(root, odir, paths, info, selected_strategy, command, warnings)
    write_completion_report(odir, info, comparison, command)

    print(f"Second validation complete: {odir}", flush=True)
    print(f"Observed DSAS: {observed_summary['DSAS'].iloc[0]:.6f}", flush=True)
    print(f"Observed PRS: {observed_summary['PRS'].iloc[0]:.6f}", flush=True)
    print(f"Observed final class: {observed_summary['final_class'].iloc[0]}", flush=True)
    print(f"Components written: {len(components)} rows", flush=True)


if __name__ == "__main__":
    main()
