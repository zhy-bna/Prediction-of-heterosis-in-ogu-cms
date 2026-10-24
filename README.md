# Prediction of Heterosis in Ogu-CMS Rapeseed

This repository contains the R scripts and data processing pipelines for the manuscript:
**"Integrative Genomic and Transcriptomic Analysis Uncovers Predictors of Biomass Heterosis in Rapeseed (Brassica napus L.)"**

## 1. System Requirements
*   **R** (version 4.5.2 or higher)
*   Required R packages: `data.table`, `vcfR`, `glmnet`, `rrBLUP`, `BGLR`, `kernlab`, `caret`, `pls`, `vroom`, `parallel`, `magrittr`.
*(Install using `install.packages(c("data.table", "glmnet", ...))`)*

## 2. File Structure & Scripts
The R scripts are located in the main directory of this repository:

*   **`Calculate additive and dominance PGSI per window.R`**: Calculates the additive and dominance PGSI matrices for genomic windows (100 kb, 500 kb, and 1 Mb).
*   **`PGSI_cv.R`**: Trains the PGSI-based LASSO models. It includes 10-fold CV, LOMO-CV, LORO-CV, external test set evaluation, and **crucially, the 100-iteration Bootstrap stability analysis (with replacement) to identify stable predictive windows (selection frequency ≥ 70%)**.
*   **`LASSO_GBLUP_SNP_cv.R`**: Benchmarks SNP-based models (LASSO_SNP, GBLUP, GBLUP_AD) using 5% random SNP subsampling across 50 seeds. **Note: A standard two-step approach (LASSO selection followed by OLS refitting) was applied to ensure fair comparison.**
*   **`BayesA_cv.R`**, **`BayesB_cv.R`**, **`PLS_cv.R`**, **`RKHS_cv.R`**, **`SVM_cv.R`**: Scripts for other benchmark genomic prediction models.
*   **`prediction for 49,550 combinations.R`**: Predicts heterosis values for 49,550 in silico hybrid combinations using a logit transformation constrained to (-100%, +100%).

## 3. Data Availability
*   **Sequencing Data**: Raw reads are available at NCBI (BioProject IDs: PRJNA1281546, PRJNA664250, PRJNA476657).
*   **Phenotype Data**: `phenotype_training_336.csv` and `phenotype_test_92.csv`.
*   **⚠️ Large Datasets (PGSI Matrices & Prediction Indices)**: 
    Due to GitHub's file size limitations, the full raw PGSI matrices (`100kb`, `500kb`, `1Mb` for additive and dominance) and the 991-accession prediction similarity index (`991_Restore.50_similarity_index.csv`, 860 MB) are **provided as Release Assets**. 
    👉 **Please download these large CSV files directly from the [Releases page](../../releases) or the Zenodo repository.**

## 4. How to Run
1. Clone or download this repository.
2. Download the large datasets from the `v1.0.0` Release page and place them in the same folder as the scripts.
3. Ensure all R packages are installed.
4. Run the scripts according to your analysis goal:
   *   For PGSI calculation: `Calculate additive and dominance PGSI per window.R`
   *   For PGSI models & Bootstrap: `PGSI_cv.R`
   *   For SNP benchmarks: `LASSO_GBLUP_SNP_cv.R`, `BayesA_cv.R`, etc.
   *   For 49,550 in silico predictions: `prediction for 49,550 combinations.R`
