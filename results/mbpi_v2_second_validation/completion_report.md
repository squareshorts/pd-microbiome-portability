# MBPI-v2 Second Validation Completion Report

## Selected Dataset

Selected dataset: MetaGenoPolis / Recherche Data Gouv colorectal cancer benchmark (https://doi.org/10.57745/7IVO3E).

It was selected because it is disease-independent from PD, has public multi-cohort shotgun metagenomic abundance data, explicit CRC/control labels, explicit study labels, and minimal harmonization burden.

## Rejected Candidate Datasets

- IBD 15-dataset benchmark: rejected for this task because it is methodologically relevant but slower to retrieve and harmonize, and CD/UC/IBD subtype handling would require more phenotype decisions.
- curatedMetagenomicData disease-specific pulls: rejected for this task because it is an excellent standardized source but requires Bioconductor object retrieval and disease-specific label assembly, making it slower than the direct CRC benchmark.

## Files Created

- `dataset_feasibility_report.md`
- `second_validation_matrix.csv`
- `second_validation_metadata.csv`
- `second_validation_input_validation.csv`
- `second_validation_mbpi_v2_observed_summary.csv`
- `second_validation_mbpi_v2_components.csv`
- `second_validation_mbpi_v2_audit_table.csv`
- `second_validation_mbpi_v2_bootstrap_replicates.csv`
- `second_validation_mbpi_v2_simulation_results.csv`
- `second_validation_mbpi_v2_simulation_summary.csv`
- `second_validation_mbpi_v2_ablation_summary.csv`
- `second_validation_mbpi_v2_run_log.txt`
- `fig_second_validation_mbpi_v2_components.png`
- `fig_second_validation_mbpi_v2_components.pdf`
- `fig_second_validation_mbpi_v2_bootstrap.png`
- `fig_second_validation_mbpi_v2_bootstrap.pdf`
- `fig_second_validation_mbpi_v2_simulation_calibration.png`
- `fig_second_validation_mbpi_v2_simulation_calibration.pdf`
- `fig_second_validation_mbpi_v2_ablation.png`
- `fig_second_validation_mbpi_v2_ablation.pdf`
- `fig_second_validation_pooled_vs_loco.png`
- `fig_second_validation_pooled_vs_loco.pdf`
- `mbpi_v2_cross_dataset_comparison.csv`

## Second Validation MBPI-v2 Results

- Samples: 1850 (880 CRC, 970 controls).
- Cohorts: 12.
- Features: 413 CLR-transformed species features.
- DSAS: 0.8348 (0.8348 to 0.8348).
- PRS: 1.0000 (1.0000 to 1.0000).
- Equal-weight PRS: 0.4502.
- Final class: non-portable.
- DSAS calibration AUROC: 1.0000.
- PRS calibration AUROC: 0.7700.
- Gated MBPI-v2 AUROC: 0.8720.
- Null false non-portability/cohort-conditioned rate: 0.0000.

## Comparison With PD Validation Case

- PD validation: DSAS 0.8625, PRS 1.0000, final class `non-portable`.
- CRC second validation: DSAS 0.8348, PRS 1.0000, final class `non-portable`.

The second dataset strengthens the methodological claim that MBPI-v2 can be applied beyond Parkinson's disease because the full pipeline completed on an independent non-PD disease area with explicit cohort labels and binary disease/control labels. It should still be described as one additional successful validation case, not proof of universal generality.

## Exact Command To Reproduce

```bash
python scripts/mbpi_v2/run_second_validation.py --root . --simulation-reps 5 --bootstrap-B 50 --seed 20260610 --prevalence-threshold 0.20
```

## Missing Inputs Or Limitations

- No manuscript text was modified.
- Adenoma samples were excluded rather than merged with CRC or controls.
- Cohorts with one class or fewer than 10 samples per class were excluded from the binary LOCO validation.
- The selected feature level is METEOR MSP species identifiers rather than genus labels.
- Bootstrap and simulation replicate counts are recorded in `second_validation_mbpi_v2_run_log.txt`; increase them for publication-grade interval stability if runtime allows.
- In this scoped run, selected PRS AUROC was 0.7700 and heterogeneity-alone AUROC was 0.7917; interpret the PRS calibration as a completed second-dataset application, not as evidence that calibrated PRS improved over every ablation.
- Bootstrap intervals for DSAS and PRS were degenerate in this run because the bootstrapped disease-signal and calibrated PRS values were saturated at the observed values.
