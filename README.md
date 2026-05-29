# Cross-Species RNA-seq Biomarker Pipeline

![R](https://img.shields.io/badge/R-4.3%2B-blue)
![Platform](https://img.shields.io/badge/Linux-WSL-green)
![Pipeline](https://img.shields.io/badge/Workflow-Automated-orange)
![Validation](https://img.shields.io/badge/Validation-External-success)

---
## Project Overview

This repository contains a modular RNA-seq analysis workflow for **cross-species biomarker discovery**.  
The project uses a mouse transcriptomic discovery cohort to identify candidate genes, transfers them to human through **strict offline ortholog mapping**, and validates them in human liver biopsy cohorts.

The workflow is designed not only for analysis, but also for **engineering-style delivery**:
- config-driven execution
- structured outputs
- clear logging
- fault tolerance for unstable small-sample modeling
- independent external validation

This project is suitable for demonstrating competencies required in **bioinformatics analyst / engineer / R&D** roles, especially for positions involving **RNA-seq, biomarker discovery, NGS workflow development, and cross-dataset validation**.

---
## My Contributions

I independently completed the following work in this project:

- Designed and implemented the full cross-species RNA-seq biomarker workflow in **R under Linux/WSL**
- Wrote modular scripts for:
  - differential expression analysis
  - functional enrichment
  - mouse-side candidate feature screening
  - human cohort validation
- Built a **strict offline ortholog mapping strategy** to replace unstable fallback conversion logic
- Refactored the human validation stage from an internal train/test split into a **development + external validation** design
- Standardized the project using:
  - YAML configuration
  - fixed directory layout
  - consistent output naming
  - explicit logging
- Summarized the project into a recruiter-friendly and technically reviewable format

---

## Key Features

### 1. End-to-End Automation
The pipeline runs from raw count inputs to final human validation outputs with modular R scripts and a shell entry point.

### 2. Strict 1:1 Ortholog Mapping
Mouse candidate genes are transferred to human using a **strict ortholog reference table**, avoiding loose symbol conversion or `toupper()` fallback logic.

### 3. Fault-Tolerant Workflow
When modeling conditions are unstable, the workflow preserves intermediate outputs and supports downgrade to candidate-mode processing instead of failing silently.

### 4. Reproducibility
The project uses:
- central YAML config
- deterministic output folders
- explicit file naming
- saved feature scaling parameters
- session info logs

---

## Directory Structure

```text
project_root/
├── config/
│   ├── config.yaml
│   └── sample_info.tsv
├── data/
│   └── reference/
│       └── mouse_human_ortholog_strict.tsv
├── logs/
├── results/
│   ├── 01_DESeq2/
│   ├── 02_Enrichment/
│   └── 03_MachineLearning/
│       ├── mouse_feature_selection/
│       └── human_validation/
├── scripts/
│   ├── utils/
│   │   ├── utils_build_ortholog_map.R
│   │   └── utils_convert_mouse_human_ids.R
│   ├── 00_fetch_and_qc.sh
│   ├── 00_upstream_align_quant.sh
│   ├── 01_deseq2_analysis.R
│   ├── 02_enrichment.R
│   ├── 03_lasso_mouse.R
│   ├── 04_1_hv_common.R
│   ├── 04_2_hv_mapping_io.R
│   ├── 04_3_hv_cohort_io.R
│   ├── 04_4_hv_modeling.R
│   ├── 04_5_hv_calibration_threshold.R
│   ├── 04_6_hv_recalibration.R
│   ├── 04_human_validation.R
│   ├── check_inputs.R
│   └── check_pipeline_health.py
├── setup/
├── .gitignore
├── environment_rnaseq_core.yml
├── README.md
└── run_pipeline.sh
```