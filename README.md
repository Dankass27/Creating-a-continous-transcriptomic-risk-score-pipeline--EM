# Development and Validation of a Continuous Transcriptomic Risk Score in Paediatric Medulloblastoma

## Overview

This repository contains the R code developed for my **BSc Biomedical Science dissertation project**, completed between **September 2025 and May 2026**.

The project investigated whether transcriptomic data could be used to develop **continuous prognostic risk scores for paediatric medulloblastoma**, with separate models developed for:

* **Group 3 / Group 4 medulloblastoma**
* **SHH medulloblastoma**

Rather than assigning patients only to predefined molecular or clinical risk groups, the aim was to generate a continuous score representing an individual patient's predicted risk based on tumour gene-expression patterns.

The project combined **transcriptomic analysis, survival statistics, penalised regression and pathway enrichment analysis**, with additional validation performed using an independent RNA-sequencing dataset.

---

## Project Aim

The primary aim was to:

> Develop and validate continuous transcriptomic risk-score models capable of stratifying patients with paediatric medulloblastoma according to overall survival.

Separate modelling approaches were used for Group 3/4 and SHH tumours because of the substantial biological and prognostic differences between medulloblastoma molecular subgroups.

---

## Analysis Workflow

The analysis was conducted primarily in **R** and consisted of several stages.

### 1. Data Preparation

Publicly available transcriptomic and clinical data were processed and combined to produce analysis-ready datasets containing:

* Gene-expression measurements
* Molecular subgroup information
* Patient survival data
* Clinical characteristics
* Histological information
* Metastatic status
* MYC/MYCN amplification status where available

The primary discovery cohort was derived from **GSE85217**.

---

### 2. Exploratory Transcriptomic Analysis

Dimensionality-reduction approaches were used to investigate transcriptomic structure within the dataset.

These included:

* **Principal Component Analysis (PCA)**
* **UMAP**

These analyses were used to visualise relationships between molecular subgroups and explore broader patterns within the transcriptomic data before prognostic modelling.

---

### 3. Survival Analysis

Clinical and molecular variables were investigated for associations with overall survival using:

* **Kaplan–Meier survival curves**
* **Log-rank testing**
* **Univariable Cox proportional hazards regression**

Variables investigated included molecular subgroup, histology, metastatic status, age and MYC/MYCN amplification.

These analyses established the clinical and molecular survival patterns present within the discovery cohort before transcriptomic risk modelling.

---

## Transcriptomic Risk-Score Development

Penalised Cox proportional hazards regression was used to develop transcriptomic prognostic models.

Three regularisation approaches were explored:

* **Ridge regression**
* **Elastic net regression**
* **LASSO regression**

Model optimisation was performed using **10-fold cross-validation**.

Separate models were developed for **Group 3/4** and **SHH** medulloblastoma.

### Group 3 / Group 4 Model

The final Group 3/4 elastic-net model contained **22 genes** and achieved a discovery-cohort concordance index of approximately **0.77**.

The resulting continuous risk score demonstrated a biological gradient between lower-risk, predominantly Group 4-like tumours and higher-risk, predominantly Group 3-like tumours.

High-risk tumours were associated with recognised adverse biological characteristics, including enrichment for features such as **MYC amplification and large-cell/anaplastic histology**.

### SHH Model

The final SHH elastic-net model contained **30 genes** and demonstrated strong prognostic discrimination within the discovery cohort, with a concordance index of approximately **0.93**.

---

## Risk Stratification

Although the models generated a **continuous risk score**, patients were additionally divided into risk groups using score tertiles to aid interpretation and visualisation.

Patients were classified as:

* Low risk
* Intermediate risk
* High risk

Kaplan–Meier analysis was then used to assess differences in overall survival between these transcriptomically derived risk groups.

---

## External Validation

To assess whether the identified prognostic signal could generalise beyond the original microarray dataset, the Group 3/4 model was projected onto an independent **RNA-sequencing validation cohort**.

This represented an important test of model robustness because the discovery and validation datasets were generated using different transcriptomic platforms.

The externally validated Group 3/4 model retained prognostic discrimination, with a corrected concordance index of approximately **0.64**, and transcriptomic risk tertiles remained significantly associated with survival.

---

## Functional Gene Set Enrichment Analysis

**Fast Gene Set Enrichment Analysis (FGSEA)** was performed to investigate the molecular biology associated with different regions of the continuous risk spectrum.

This analysis was used to identify biological pathways and transcriptional programmes associated with higher- and lower-risk tumour profiles.

The results supported the interpretation of medulloblastoma risk as a **biological continuum**, particularly within Group 3/4 disease, rather than simply a set of discrete prognostic categories.

---

## Repository Contents

The repository contains both the original analytical workflow and subsequently refined versions of the code.

### Source Code

Contains the full analysis workflow developed throughout the dissertation project, including exploratory analyses, intermediate modelling steps and supporting analyses.

### Condensed Source Code

A refined and organised version of the main analytical pipeline.

This contains the principal stages required for:

* Data preparation
* Exploratory transcriptomic analysis
* Survival analysis
* Penalised Cox modelling
* Risk-score generation
* Risk stratification
* Model evaluation
* Cross-dataset validation

### FGSEA Analysis

Contains the code used for functional gene set enrichment analysis and biological interpretation of the transcriptomic risk models.

### Tables

Contains code used to generate statistical, clinical and model-performance tables presented within the dissertation.

### Supporting / Appendix Analyses

Additional scripts contain supporting analyses and figures used during the dissertation and supplementary analysis process.

---

## Key Methods

The project involved the application of:

* R programming
* Transcriptomic data analysis
* Gene-expression data processing
* Kaplan–Meier survival analysis
* Log-rank testing
* Cox proportional hazards regression
* Penalised Cox regression
* Ridge regression
* LASSO regression
* Elastic-net regression
* 10-fold cross-validation
* Concordance-index model evaluation
* PCA
* UMAP
* FGSEA
* Cross-platform transcriptomic validation
* Data visualisation
* Clinical and molecular data integration

---

## Key Findings

The project demonstrated that transcriptomic information could be used to generate continuous prognostic scores within medulloblastoma molecular subgroups.

In particular:

* Separate transcriptomic risk models were successfully developed for Group 3/4 and SHH medulloblastoma.
* Penalised Cox regression enabled prognostic gene-expression signatures to be derived from high-dimensional transcriptomic data.
* The Group 3/4 score captured a biological and prognostic continuum spanning predominantly Group 4-like lower-risk tumours to Group 3-like higher-risk tumours.
* The Group 3/4 model retained prognostic information following projection from a microarray discovery dataset into an independent RNA-sequencing cohort.
* Pathway enrichment analysis provided additional biological interpretation of the molecular processes associated with the risk spectrum.

These findings support the potential value of **continuous molecular risk modelling** as a complement to conventional categorical medulloblastoma risk classification.

---

## Project Context

This work was completed as my final-year dissertation for the **BSc Biomedical Science degree at Northumbria University**.

The project was designed as an academic research project and the models contained within this repository are **not intended for clinical use**.

---

## Author

**Ewan M.**

BSc Biomedical Science

Completed: **May 2026**

---

## Disclaimer

This repository contains code produced for an undergraduate research dissertation.

The risk models presented here are exploratory research models developed using retrospective publicly available datasets. They have not undergone prospective clinical validation and should not be used to inform clinical decision-making.
