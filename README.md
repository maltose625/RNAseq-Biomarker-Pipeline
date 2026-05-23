# Cross-Species RNA-seq Biomarker Pipeline

![R](https://img.shields.io/badge/R-4.3%2B-blue)
![Platform](https://img.shields.io/badge/Linux-WSL-green)
![Pipeline](https://img.shields.io/badge/Workflow-Automated-orange)
![Validation](https://img.shields.io/badge/Validation-External-success)

> **TL;DR**  
> This project is an end-to-end, config-driven cross-species RNA-seq biomarker discovery pipeline developed with **R + Linux/WSL**.  
> It starts from **mouse RNA-seq differential expression analysis**, performs **strict 1:1 mouse-human ortholog mapping**, and validates transferred biomarker candidates in **independent human clinical cohorts** using **LASSO logistic regression**.
>
> **Key outputs:**  
> - >1900 significant DE genes in mouse  
> - 1054 ortholog-transferred human features  
> - Development cohort AUC = **0.9811**  
> - External validation cohort AUC = **0.7682**

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

## Quick Facts

| Item | Details |
|------|---------|
| Project Type | Cross-species RNA-seq biomarker discovery pipeline |
| Main Language | R |
| Execution Environment | Linux / WSL |
| Configuration | YAML-based |
| Discovery Data | Mouse RNA-seq |
| Validation Data | Human GEO cohorts |
| Core Methods | DESeq2, clusterProfiler, glmnet, ROC analysis |
| Ortholog Strategy | Strict offline 1:1 mouse-human mapping |
| Development Cohort | GSE193066 (N=106) |
| External Validation Cohort | GSE193080 (N=59) |
| Final Goal | Human biomarker transfer and validation |

---

## Skills Demonstrated

This project demonstrates the following job-relevant technical and engineering skills:

- **RNA-seq differential expression analysis** using `DESeq2`
- **Functional enrichment analysis** using `clusterProfiler`
- **Cross-species biomarker transfer** using strict mouse-human ortholog mapping
- **Machine learning-based feature selection** using LASSO logistic regression (`glmnet`)
- **Independent external validation** in a second human cohort
- **Linux / WSL command-line workflow execution**
- **Config-driven reproducible pipeline design**
- **Structured logging and maintainable result outputs**
- **Fault-tolerant workflow design** for small-sample instability

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

### 4. External Validation
The human validation stage uses:
- **Development Cohort:** GSE193066, N=106
- **External Validation Cohort:** GSE193080, N=59

### 5. Reproducibility
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
├── README.md
├── run_pipeline.sh
├── .gitignore
├── config/
│   ├── config.yaml
│   └── sample_info.tsv
├── data/
│   └── reference/
│       └── mouse_human_ortholog_strict.tsv
├── logs/
├── scripts/
│   ├── check_inputs.R
│   ├── 01_deseq2_analysis.R
│   ├── 02_enrichment.R
│   ├── 03_lasso_mouse.R
│   ├── 04_human_validation.R
│   └── utils_parse_mgi_ortholog.R
└── results/
    ├── 01_DESeq2/
    ├── 02_Enrichment/
    └── 03_MachineLearning/
        ├── mouse_feature_selection/
        └── human_validation/
