# Statistical Reporting Standards for Nature

## Core Principle

Every quantitative claim must include: **test name**, **sample size (n)**, **test statistic**, **exact P value**, and **effect size**. Nature editors and reviewers will reject manuscripts with incomplete statistical reporting.

---

## P-Value Reporting

### Format Rules
| P Value Range | Format | Example |
|---|---|---|
| P ≥ 0.001 | Exact, 2 significant figures | P = 0.023 |
| 0.0001 ≤ P < 0.001 | Exact or threshold | P = 4.7 × 10⁻⁴ |
| P < 0.0001 | Scientific notation | P = 3.2 × 10⁻⁸ |
| Very small | Threshold | P < 1 × 10⁻¹⁵ |

### What to Report Alongside P Values
```
"Expression of CCL18 was significantly higher in TB-affected cells
compared to controls (mean difference = 2.3 log₂ fold-change,
95% CI [1.8, 2.7], Wilcoxon rank-sum test, n = 12,198 TB,
n = 3,049 control, P = 3.2 × 10⁻⁸)."
```

### P-Value Language
| P Value | Acceptable Language |
|---|---|
| P < 0.05 | "significantly different" (with the test name) |
| P ≥ 0.05 | "not significantly different" (never "trend toward significance") |
| P = 0.05 exactly | "borderline" — acknowledge ambiguity |

**NEVER say**: "trending toward significance" (P = 0.07). It either is or isn't.

---

## Effect Sizes

Report effect sizes for all primary outcomes:

| Comparison Type | Effect Size Measure | Interpretation |
|---|---|---|
| Two-group continuous | Cohen's d | Small: 0.2, Medium: 0.5, Large: 0.8 |
| Two-group categorical | Odds ratio (OR) or relative risk (RR) | With 95% CI |
| Correlation | Pearson r or Spearman ρ | Small: 0.1, Medium: 0.3, Large: 0.5 |
| Classification | AUC-ROC | Good: >0.8, Excellent: >0.9 |
| Regression | R², adjusted R² | Proportion of variance explained |
| Fold change | Log₂ fold change | Directionality + magnitude |

---

## Sample Size Reporting

### What Counts as n?
- **Biological replicates**: Independent samples, patients, animals
- **NOT technical replicates**: Repeated measurements of the same sample
- **For single-cell**: n = number of cells, but also report n_patients/donors

### Where to Report
- In figure legends: "n = X cells from Y donors"
- In Methods: full experimental details
- In Results text: when first presenting the data

---

## Multiple Testing Correction

When performing many statistical tests simultaneously:

| Number of Tests | Correction Method | When to Use |
|---|---|---|
| 2–5 tests | Bonferroni | Conservative, small number of tests |
| 5–100 tests | Benjamini-Hochberg (FDR) | Genomics, multiple endpoints |
| >100 tests | FDR with q-value threshold | Genome-wide studies |

**Always state**: "P values were adjusted for multiple comparisons using the Benjamini-Hochberg method (FDR < 0.05)."

---

## Error Bars and Uncertainty

### Define Error Bars Explicitly
Every figure legend with error bars must state:
- "Data are mean ± s.d." (standard deviation — variability)
- "Data are mean ± s.e.m." (standard error of mean — precision)
- "Box plots show median (line), IQR (box), and 1.5× IQR (whiskers)"

### When to Use Which
| Measure | Use When |
|---|---|
| s.d. | Describing the spread/variability of data |
| s.e.m. | Comparing means between groups |
| 95% CI | Estimating population parameters |
| IQR | Non-normal data, box plots |

---

## Machine Learning-Specific Statistics

### Classification Metrics
Report ALL of these for classification tasks:
- **AUC-ROC** with 95% CI (DeLong method or bootstrap)
- **Accuracy**, **Precision**, **Recall**, **F1-score**
- **MCC** (Matthews Correlation Coefficient) for imbalanced data
- **Confusion matrix** (in Extended Data or Supplementary)

### Cross-Validation Reporting
```
"Five-fold stratified cross-validation yielded a mean AUC-ROC of
0.947 ± 0.012 (mean ± s.d. across folds) for the stacking ensemble,
compared to 0.932 ± 0.015 for Random Forest, 0.941 ± 0.013 for
XGBoost, and 0.938 ± 0.014 for LightGBM."
```

### SHAP Values
- Report mean |SHAP| for feature importance
- Include confidence intervals from bootstrap resampling
- State the number of background samples used for interventional SHAP
- Cite Lundberg et al. (2020) Nat. Mach. Intell.

### Model Comparison
- Use DeLong test to compare AUC-ROC between models
- Report paired differences with P values
- Use McNemar test for comparing classification accuracy

---

## Reporting in Methods Section

### Template for Statistical Methods Paragraph

```
Statistical analyses were performed using Python 3.11 with
scikit-learn v1.3 (ref. X), SciPy v1.11 (ref. Y), and statsmodels
v0.14 (ref. Z). Group comparisons were performed using Wilcoxon
rank-sum tests for non-normally distributed data and two-sided
Student's t-tests for normally distributed data. Distribution
normality was assessed using the Shapiro-Wilk test. P values were
adjusted for multiple comparisons using the Benjamini-Hochberg
procedure, with a false discovery rate threshold of 0.05.
Classification performance was evaluated using five-fold stratified
cross-validation, with AUC-ROC as the primary metric. 95% confidence
intervals for AUC-ROC were computed using 1,000 bootstrap iterations.
All statistical tests were two-sided unless otherwise noted.
```

---

## Enrichment & Pathway Analysis Statistics

### Gene Set Enrichment
- Report enrichment score, normalized enrichment score (NES), P value, and FDR q-value
- State the gene set database used (MSigDB, KEGG, GO, Reactome) with version
- Report the number of gene sets tested
- For CMap/L1000: report connectivity score, enrichment score, and permutation-based P value

### Reporting Template
```
"Gene set enrichment analysis revealed significant enrichment of
inflammatory response pathways (NES = 2.3, FDR q = 0.001) and
complement activation (NES = 1.9, FDR q = 0.003) in TB-affected cells.
Analysis was performed using GSEA v4.3 against MSigDB Hallmark gene
sets (v7.5.1). FDR q-values were computed from 1,000 permutations."
```

---

## Survival Analysis Statistics (If Applicable)

### Required Elements
- **Kaplan-Meier**: Report median survival time with 95% CI for each group
- **Log-rank test**: Report chi-squared statistic and P value
- **Cox regression**: Report hazard ratio (HR) with 95% CI and P value
- **Censoring**: State the censoring mechanism and proportion censored
- **Follow-up**: Report median follow-up time

### Reporting Template
```
"Patients with high CCL18 expression (above median) had significantly
shorter overall survival (median 18.3 months, 95% CI [14.2, 22.1])
compared to low-expression patients (median 28.7 months, 95% CI
[23.5, 34.2]; log-rank test, P = 0.003). In multivariable Cox
regression adjusting for age, sex, and stage, CCL18 expression
remained independently prognostic (HR = 1.82, 95% CI [1.23, 2.71],
P = 0.003)."
```

---

## Dimensionality Reduction Statistics

### UMAP/t-SNE
- Report all non-default parameters (n_neighbors, min_dist, perplexity)
- State the number of principal components used as input
- Report the number of highly variable genes selected
- State the random seed for reproducibility

### Clustering
- Report the algorithm (Leiden, Louvain) and resolution parameter
- Report the number of clusters identified
- State the metric used (silhouette score, modularity) if applicable

---

## Nature Reporting Summary

Nature requires completion of their standardized Reporting Summary form covering:

1. **Statistics**: Tests used, multiple testing corrections, sample sizes, power analysis
2. **Study design**: Randomization, blinding, replication
3. **Data collection**: Methods, instruments, software
4. **Data analysis**: Code availability, software versions, custom algorithms
5. **Code/Data**: Repository links, accession numbers

Download: https://www.nature.com/documents/nr-reporting-summary.pdf

