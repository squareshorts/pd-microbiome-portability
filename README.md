Cross-Site Portability of Parkinson’s Disease Microbiome Classifiers

This repository contains the complete reproducible analysis pipeline supporting the manuscript:

“Site-Level Distribution Shift Dominates Disease Signal and Undermines Cross-Site Generalization in High-Dimensional Microbiome Classifiers.”

Overview

This project evaluates the stability of microbiome-based Parkinson’s disease (PD) classifiers under cross-site distribution shift across three geographically distinct cohorts (Finland, Malaysia, USA; total n = 682).

The study demonstrates that:

Within-cohort cross-validation yields high discrimination (AUC up to 0.944).

Leave-one-cohort-out validation produces substantial degradation (AUC as low as 0.506).

Cohort identity explains 68.9% of compositional variance.

PD status explains 0.92%.

Disease-associated effect vectors are nearly orthogonal to the dominant ecological axis (PC1).

These results formalize cross-site portability as a structural limitation in high-dimensional compositional biomedical modeling.

Repository Structure
data/processed/ipd_paper2/
    Harmonized genus-level matrices and metadata

scripts/R/paper2/
    Full reproducible analysis pipeline (C21–C33)

results/paper2/
    Summary tables, fold-level outputs, LOCO predictions, figures
Reproducibility

R version: 4.4.2

Deterministic random seeds fixed

All performance tables generated from scripts

All figures reproducible from included code

To reproduce the full analysis:

source("scripts/R/paper2/C21_clr_shared2.R")
...
source("scripts/R/paper2/C33_pd_effect_alignment_angles.R")

Scripts are numerically ordered to match manuscript sections.

License

Specify your chosen license here (e.g., MIT, BSD-3-Clause, GPL-3.0).

DOI

Zenodo archive:
https://doi.org/10.5281/zenodo.18674302