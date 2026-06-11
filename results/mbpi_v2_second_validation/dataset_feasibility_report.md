# MBPI-v2 Second Validation Dataset Feasibility Report

## Decision

Selected dataset: MetaGenoPolis / Recherche Data Gouv colorectal cancer benchmark.

Reason: it is the fastest credible independent validation case because it provides a single public repository with manually curated metadata, explicit study labels, explicit CRC/healthy/adenoma phenotypes, and ready-to-use species abundance tables. It also has strong methodological value for MBPI-v2 because the retained CRC/control subset contains 1850 samples from 12 independent studies after transparent quality and cohort-balance filters.

## Candidate 1: MetaGenoPolis / Recherche Data Gouv CRC benchmark

- Dataset name: Taxonomic profiles, functional profiles and manually curated metadata of human fecal metagenomes from public projects coming from colorectal cancer studies.
- Disease area: Colorectal cancer.
- Source/accession or repository: https://doi.org/10.57745/7IVO3E; selected files https://doi.org/10.57745/DQBQD3 and https://doi.org/10.57745/LCAR4M.
- Number of cohorts/studies: 15 available; 12 retained for binary CRC/control LOCO validation.
- Approximate sample size: 2,340 total metadata rows; 1850 retained after excluding adenoma, quality-flagged rows, and cohorts without at least 10 CRC and 10 controls.
- Feature type: METEOR metagenomic species pangenome abundance, CLR transformed after a 0.20 non-zero prevalence filter (413 retained features).
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
