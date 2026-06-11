# MBPI-v2 Algorithm Specification

MBPI-v2 is a two-stage diagnostic algorithm for microbiome biomarker portability. It separates disease-signal adequacy from portability risk. This is a technical specification, not manuscript text.

## Inputs

- CLR-transformed microbiome matrix `X` with dimension `n x p`
- phenotype vector `y`
- cohort vector `c`
- optional classifier summaries
- optional genus-level FDR and fixed-effect meta-analysis tables

Observed run primary input: `C:\work\pd-microbiome-portability\results\portability_analysis\dataset_clr.csv`

## Stage 1: Disease Signal Adequacy

The Disease Signal Adequacy Score, DSAS, is bounded in `[0, 1]`.

Components:

- `DSAS_pooled_auroc = sqrt(max(0, median(AUROC_pooled) - 0.5) / 0.5)`
- `DSAS_loco_auroc = sqrt(max(0, median(AUROC_LOCO) - 0.5) / 0.5)`
- `DSAS_pd_main_fdr = min(1, proportion(PD-main BH q < 0.05) / 0.25)`
- `DSAS_fixed_effect_fdr = min(1, proportion(fixed-effect meta-analysis BH q < 0.05) / 0.20)`
- `DSAS_effect_size = min(1, median absolute PD beta among PD-main q < 0.05 / 0.25)`

Weights:

DSAS_pooled_auroc        0.50
DSAS_loco_auroc          0.10
DSAS_pd_main_fdr         0.20
DSAS_fixed_effect_fdr    0.15
DSAS_effect_size         0.05

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

Equal-weight PRS and calibrated-weight PRS are both reported. For logistic calibration strategies, PRS is the calibrated logistic probability from non-negative component coefficients; normalized component weights are also reported for interpretability. For manual and equal-weight strategies, PRS is a bounded weighted mean. The primary calibrated strategy selected for this run is `nonnegative_logistic`.

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
python scripts/mbpi_v2/run_mbpi_v2_pipeline.py --simulation-reps 50 --bootstrap-B 500 --seed 20260608
```
