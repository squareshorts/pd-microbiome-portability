# Biomarker Portability Analysis Technical Report

## 1. Dataset Overview

The analysis uses the harmonised genus-level dataset representing the true manuscript cohort.

| Metric | Value |
| --- | --- |
| Total Samples | 683 |
| Number of Genus Features (CLR) | 62 |
| Cohort: Finland - Total | 148 |
| Cohort: Finland - PD | 74 |
| Cohort: Finland - Control | 74 |
| Cohort: Malaysia - Total | 198 |
| Cohort: Malaysia - PD | 101 |
| Cohort: Malaysia - Control | 97 |
| Cohort: USA - Total | 337 |
| Cohort: USA - PD | 204 |
| Cohort: USA - Control | 133 |

## 2. Disease-Axis Geometry & Variance

- **PC1 Variance Explained:** 61.2%
- **Cosine Similarity (PC1 vs PD Vector):** 0.067
- **Angle:** 86.1°
- **CDDR (Cohort/PD Ratio from PERMANOVA):** 258.7
  - Cohort R2: 0.7505
  - PD R2: 0.0029

## 3. Batch Correction & Domain Adaptation

**DEBIAS-M Evaluation**
DEBIAS-M is not valid for strict LOCO under the available workflow without held-out phenotype-label leakage.

**Alternative: Transductive Cohort-Centering**
Transductive cohort-centering failed to produce a valid generalization of the PD signature. Removing the cohort effect transductively resulted in degenerate classifiers with AUROC = 0.500, sensitivity = 1.0, specificity = 0.0, indicating that any successful pooled prediction was leveraging cohort bias rather than true disease biology.

## 4. Model Performance Summary (Baseline)

| Model | Pooled AUROC | Mean LOCO AUROC | BPG AUROC | Pooled AUPRC | Mean LOCO AUPRC | BPG AUPRC |
| --- | --- | --- | --- | --- | --- | --- |
| ElasticNet | 0.810 | 0.600 | 0.210 | 0.821 | 0.648 | 0.173 |
| RandomForest | 0.881 | 0.515 | 0.366 | 0.894 | 0.573 | 0.321 |
| XGBoost | 0.919 | 0.547 | 0.371 | 0.926 | 0.601 | 0.325 |

## 5. Model Performance Summary (Transductive Corrected)

| Model | Pooled AUROC | Mean LOCO AUROC | Mean LOCO Sens | Mean LOCO Spec |
| --- | --- | --- | --- | --- |
| ElasticNet | 0.500 | 0.500 | 1.000 | 0.000 |
| RandomForest | 0.500 | 0.500 | 1.000 | 0.000 |
| XGBoost | 0.500 | 0.500 | 1.000 | 0.000 |

## 6. Feature Stability Summary

| Model | Mean Overlap Top 10 | Mean Overlap Top 20 | Mean Jaccard Top 20 |
| --- | --- | --- | --- |
| ElasticNet | 2.3 | 10.0 | 0.334 |
| RandomForest | 3.0 | 11.0 | 0.380 |
| XGBoost | 2.0 | 7.3 | 0.225 |

## 7. Generated Figures

- **Figure ML1:** `Figure_ML1_Pooled_vs_LOCO.png` (Pooled CV vs LOCO metrics)
- **Figure ML2:** `Figure_ML2_LOCO_by_Cohort.png` (LOCO performance broken down by held-out cohort)
- **Figure ML3:** `Figure_ML3_BPG.png` (Biomarker Portability Gap: Pooled - LOCO)
- **Figure ML4:** `Figure_ML4_PCA.png` (PCA of Raw vs Global Residualized CLR)
- **Figure ML5:** `Figure_ML5_CDDR_Geometry.png` (CDDR and Geometry summary)
- **Figure ML6:** `Figure_ML6_Feature_Stability.png` (Top-10 and Top-20 feature stability across cohorts)
