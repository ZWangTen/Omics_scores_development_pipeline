# 📊 Multi-Omics Elastic Net Pipeline (R + targets)

A reproducible and scalable R pipeline for **high-dimensional omics data preprocessing, feature selection, and risk score construction using elastic net. Built with the [`targets`](https://docs.ropensci.org/targets/) framework for workflow management, this pipeline supports **data quality control, normalization, imputation, feature filtering, and model selection** in a fully reproducible manner.

## 🚀 Overview

This pipeline is designed for analyzing high-dimensional biological data (e.g., metabolomics, proteomics, methylation), integrating phenotype data, and identifying predictive features associated with clinical outcomes.

### Key features:

- Reproducible workflow using `targets`
- Data QC: normality + missingness assessment  
- Optional preprocessing: rank-normalization + imputation  
- Feature filtering based on association testing  
- Elastic net model with hyperparameter tuning  
- Modular, scalable, and easy to extend  

## 🧬 Workflow

1. Data Input  
2. Quality Control  
3. Preprocessing (optional)  
4. Feature Filtering (optional)  
5. Modeling (Elastic Net)  
6. Output generation  

## ⚙️ Installation

```bash
git clone https://github.com/ZWangTen/Omics_scores_development_pipeline.git
cd YOUR_REPO
```

```r
install.packages(c("targets", "tidyverse", "glmnet"))
```

## ▶️ Usage

```r
library(targets)
tar_make()
```

- Check `Example_run_pipeline.R` for example codes running the pipeline.
- Check `Helper_Functions_Documentation.docx` for function and arguments usage.
 
## 🔧 Configuration

Modify parameters in `_targets.R`:

- filtering (TRUE/FALSE)
- run_RN (normalization)
- run_IMP (imputation)
- alpha values

## 📥 Input

- data.csv: feature matrix with `id`
- pheno.csv: phenotype + covariates

## 📤 Output

- check_results.csv
- Selected_features.csv

## 🧠 Methods

- Rank normalization
- Half-minimum imputation
- Elastic net (glmnet)

## 👤 Author

Ziqing (Leslie) Wang  
GitHub: https://github.com/ZWangTen

## 📄 License

MIT License
