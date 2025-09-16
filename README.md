
# Breast Cancer Prediction using Bayesian Logistic Regression & GAMs

Master of Analytics project exploring **Bayesian Logistic Regression** and **Generalised Additive Models (GAMs)** for predicting breast cancer diagnosis from **Fine Needle Aspiration (FNA)**-derived features.

- **Authors:** Ali Rajabimehr (with teammate **Craig Manning**)
- **Supervisors:** **Victor Miranda**, **Patricio Andres Maturana Russel**
- **Institution:** Auckland University of Technology (AUT)

## 🔍 Overview
- **Objective:** Identify significant predictors of malignancy and compare Bayesian vs GAM performance.  
- **Dataset:** Breast Cancer Wisconsin (Diagnostic) – 569 cases, 30 features.  
- **Tools:** R (`rjags`, `coda`, `mgcv`, `car`, `corrplot`, `ggplot2`, `dplyr`, `readr`).  
- **Headline results:** GAM ≈ **94.9%** accuracy; Bayesian with informative priors ≈ **93.1%**.

## 📂 Repository Structure
```
breast-cancer-prediction-bayes-gam/
├── README.md
├── .gitignore
├── scripts/
│   └── breast_cancer_prediction.R
├── data/
│   └── data_description.md          # Put raw CSV as data/data.csv (not tracked)
├── results/
│   ├── .gitkeep
│   └── figures/
│       └── .gitkeep
└── report/
    ├── Presentation.pdf
    └── Full_Report.pdf
```

## ▶️ How to Run (R)
1. Install R packages if needed:
   ```r
   install.packages(c("dplyr","readr","ggplot2","corrplot","car","coda","mgcv","rjags","kableExtra"))
   ```
2. Put the dataset CSV at `data/data.csv` (see `data/data_description.md`).
3. From the project root, run:
   ```r
   source("scripts/breast_cancer_prediction.R")
   ```
4. Outputs will be saved to `results/` (metrics CSVs, convergence logs, figures).

## 🧠 Methods
- **Feature selection:** Correlation + ANOVA + iterative VIF pruning.
- **Bayesian Logistic Regression:** JAGS with multiple prior strengths; best model selected via **BIC**; posterior inference and **HPD** credible intervals.
- **GAM:** `mgcv` with smoothing terms per selected predictor to model non-linear relationships.

## 📊 Expected Outputs
- `results/model_comparison.csv` (Bayesian vs GAM metrics)
- `results/bayesian_metrics.csv`, `results/gam_metrics.csv`
- `results/bayes_convergence.txt`
- `results/figures/gam_smooths.png`, `results/figures/correlation_heatmap.png`

## 📢 Acknowledgements
This work was completed as part of the **Master of Analytics** at AUT, supervised by **Victor Miranda** and **Patricio Andres Maturana Russel**.  
Team collaboration with **Craig Manning**.

---

If you use or extend this repo, please star it and feel free to open issues/PRs.
