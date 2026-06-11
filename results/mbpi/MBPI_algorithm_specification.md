# Microbiome Biomarker Portability Index (MBPI)

This technical note defines the Microbiome Biomarker Portability Index (MBPI) as a reproducible diagnostic algorithm for microbiome biomarker non-portability. It is an algorithm specification, not manuscript prose.

## Inputs

- `X`: CLR-transformed microbiome matrix with dimension `n x p`.
- `y`: binary phenotype vector, where `1` denotes PD and `0` denotes control.
- `c`: cohort vector.
- Optional classifier outputs for pooled and leave-one-cohort-out (LOCO) AUROC.
- Optional feature-importance stability outputs.
- Optional genus-level p-value/FDR tables.

The observed PD validation run uses `C:\work\pd-microbiome-portability\results\portability_analysis\dataset_clr.csv` as the primary CLR and metadata source.

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
python scripts/mbpi/run_mbpi_pipeline.py --bootstrap-B 500 --bootstrap-mode fast --simulation-reps 50 --seed 20260608
```

The run log records input files, sample counts, component values, bootstrap intervals, simulation settings, ablation summaries, warnings, and file hashes where feasible.
