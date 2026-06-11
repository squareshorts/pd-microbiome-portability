# Cross-Site Portability of Parkinson's Disease Microbiome Classifiers

This repository contains the complete reproducible analysis pipeline supporting the manuscript:

**Site-Level Distribution Shift Dominates Disease Signal and Undermines Cross-Site Generalization in High-Dimensional Microbiome Classifiers.**

## Overview

This project evaluates the stability of microbiome-based Parkinson's disease (PD) classifiers under cross-site distribution shift across three geographically distinct cohorts (Finland, Malaysia, USA; total n = 682).

The study demonstrates that:

- Within-cohort cross-validation yields high discrimination (AUC up to 0.944).
- Leave-one-cohort-out validation produces substantial degradation (AUC as low as 0.506).
- Cohort identity explains 68.9% of compositional variance.
- PD status explains 0.92%.
- Disease-associated effect vectors are nearly orthogonal to the dominant ecological axis (PC1).

These results formalize cross-site portability as a structural limitation in high-dimensional compositional biomedical modeling.

## Repository Structure

```text
data/processed/ipd_paper2/
    Harmonized genus-level matrices and metadata

scripts/R/paper2/
    Main R analysis pipeline (C21-C33)

scripts/R/
    Supporting R preprocessing, sensitivity, and figure scripts

scripts/mbpi/ and scripts/mbpi_v2/
    MBPI analysis utilities and validation workflows

results/paper2/ and results/paper2_portability*/
    Main output tables, predictions, calibration outputs, and figures

manuscript/
    LaTeX manuscript source
```

## Reproducibility

- R version: 4.4.2
- Deterministic random seeds are fixed in the analysis scripts.
- Performance tables are generated from scripts.
- Figures are reproducible from the included code and processed data.

To reproduce the main paper 2 analysis, run the numbered scripts in order:

```r
source("scripts/R/paper2/C21_clr_shared2.R")
source("scripts/R/paper2/C22_ecology_structure.R")
source("scripts/R/paper2/C23_pca_plot.R")
source("scripts/R/paper2/C24_loco_logreg.R")
source("scripts/R/paper2/C25_loco_rf.R")
source("scripts/R/paper2/C26_within_cohort_logreg.R")
source("scripts/R/paper2/C27_within_cohort_rf.R")
source("scripts/R/paper2/C28_portability_gap.R")
source("scripts/R/paper2/C29_loco_threshold_sweep.R")
source("scripts/R/paper2/C30_pc1_loadings.R")
source("scripts/R/paper2/C31_pd_alignment_pc1.R")
source("scripts/R/paper2/C32_structural_asymmetry_summary.R")
source("scripts/R/paper2/C33_pd_effect_alignment_angles.R")
```

Scripts are numerically ordered to match manuscript analyses.

## Data Availability

All sequencing datasets analysed in this study are publicly available from the European Nucleotide Archive (ENA) and the NCBI Sequence Read Archive (SRA) under the following BioProject identifiers: PRJNA494620 (Malaysia), PRJEB4927 (Finland), and PRJEB14674 (United States). Corresponding study identifiers are SRP199959, ERP004264, and ERP016332, respectively.

No new sequencing data were generated. The harmonized genus-level abundance matrix, processed data tables, and all R scripts required to reproduce the analyses and figures are publicly available at Zenodo:

https://doi.org/10.5281/zenodo.18644851

## License

No license has been specified in this repository.
