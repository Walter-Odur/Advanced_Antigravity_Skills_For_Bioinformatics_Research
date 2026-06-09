---
name: scrna-drug-repurposing
description: End-to-end single-cell RNA-seq ensemble ML pipeline for disease signature extraction and computational drug repurposing. Covers scRNA-seq QC (Scanpy), ensemble classification (RF, XGBoost, LightGBM, Stacking), SHAP interpretability, CMap/L1000 drug repurposing, and publication-quality figure generation.
metadata:
    skill-author: Walter Odur
    version: "1.0.0"
    domain: bioinformatics, drug-discovery, machine-learning
risk: low
source: custom
---

# scRNA-seq Ensemble ML Drug Repurposing Pipeline

## Overview

This skill provides comprehensive guidance for building an **end-to-end computational drug repurposing pipeline** from single-cell RNA sequencing (scRNA-seq) data using ensemble machine learning. The pipeline identifies disease-driving genes via SHAP analysis and queries the Connectivity Map (CMap/L1000) to find FDA-approved drugs that may reverse the disease signature.

All methodological choices are grounded in high-impact peer-reviewed literature from *Nature*, *Cell*, *Nature Methods*, *Genome Biology*, *Nature Machine Intelligence*, and *Nature Reviews Drug Discovery*. Each pipeline step offers both a standard approach and advanced options for maximum rigor.

The canonical reference implementation is a TB granuloma drug repurposing project (see `examples/tb_granuloma/`).

## When to Use This Skill

Use this skill when:

- Building a drug repurposing pipeline from scRNA-seq data
- Performing binary classification (disease vs. control) on single-cell expression data
- Training ensemble ML classifiers (Random Forest, XGBoost, LightGBM, Stacking)
- Extracting SHAP-derived disease gene signatures from tree-based models
- Querying CMap/L1000 for candidate host-directed therapeutics
- Generating publication-quality multi-panel manuscript figures
- Working with Broad Institute Single Cell Portal (SCP) or GEO scRNA-seq datasets

## Do Not Use This Skill When

- The data is bulk RNA-seq (use DESeq2/edgeR differential expression instead)
- The task is unsupervised clustering or cell type annotation only (use the `scanpy` skill)
- The goal is deep learning on images or text (use the `ml-engineer` skill)
- You need general scikit-learn guidance without the scRNA-seq context (use `scikit-learn` skill)

## Pipeline Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                    INPUT: scRNA-seq Data                     │
│         (count matrix + cell metadata from SCP/GEO)         │
└──────────────────────┬──────────────────────────────────────┘
                       │
                       ▼
┌──────────────────────────────────────────────────────────────┐
│  STEP 1: Load & QC  (01_load_and_qc.py)                     │
│  • Load count matrix + metadata (handle SCP format quirks)   │
│  • Orient matrix (cells × genes)                             │
│  • [Advanced] Ambient RNA removal: CellBender / SoupX        │
│  • [Advanced] Doublet detection: scDblFinder / Scrublet      │
│  • Filter: min_genes=200, min_cells=3, MT% < 20             │
│  • Normalize (10K target or scran/SCTransform), select HVGs  │
│  • [Advanced] Batch correction: Harmony / scVI               │
│  OUTPUT: processed AnnData (.h5ad)                           │
└──────────────────────┬───────────────────────────────────────┘
                       │
                       ▼
┌──────────────────────────────────────────────────────────────┐
│  STEP 2: Prepare ML Data  (02_prepare_ml_data.py)            │
│  • Map disease groups → binary labels (1=disease, 0=control) │
│  • Exclude confounding groups (e.g., HIV-only)               │
│  • Extract dense feature matrix from HVG-filtered AnnData    │
│  • Stratified 80/20 train/test split                         │
│  OUTPUT: train_data.npz, test_data.npz, gene_names.csv       │
└──────────────────────┬───────────────────────────────────────┘
                       │
                       ▼
┌──────────────────────────────────────────────────────────────┐
│  STEP 3: Ensemble ML Training  (03_ensemble_ml.py)           │
│  • Train 4 classifiers:                                      │
│    1. Random Forest (Breiman 2001)                            │
│    2. XGBoost (Chen & Guestrin, KDD 2016)                     │
│    3. LightGBM (Ke et al., NeurIPS 2017)                      │
│    4. Stacking (Wolpert, Neural Networks 1992)                │
│  • 5-fold stratified CV / [Advanced] nested CV                │
│  • [Advanced] Bayesian optimization via Optuna (TPE sampler)  │
│  • [Advanced] Probability calibration (Platt / Isotonic)      │
│  • DeLong test / McNemar's test for model comparison          │
│  • Metrics: AUC-ROC, PR-AUC, MCC, Brier, F1, Cohen's κ       │
│  OUTPUT: *.joblib models, cv_results.csv, test_results.csv    │
└──────────────────────┬───────────────────────────────────────┘
                       │
                       ▼
┌──────────────────────────────────────────────────────────────┐
│  STEP 4: Model Improvement  (04_model_improvement.py)        │
│  • Expedited hyperparameter tuning (pre-selected grids)      │
│  • Learning curves for best model                            │
│  • Baseline vs. tuned comparison                             │
│  OUTPUT: *_tuned.joblib, tuned_vs_baseline_comparison.csv     │
└──────────────────────┬───────────────────────────────────────┘
                       │
                       ▼
┌──────────────────────────────────────────────────────────────┐
│  STEP 5: SHAP Analysis  (05_shap_analysis.py)                │
│  • TreeExplainer: interventional mode (Lundberg, NMI 2020)   │
│  • Extract class-1 SHAP values (disease class)               │
│  • [Advanced] SHAP interaction values for gene co-regulation  │
│  • Rank genes by mean |SHAP| → disease signature             │
│  • [Advanced] Bootstrap stability assessment (100 resamples)  │
│  • [Advanced] Pathway enrichment validation (GSEA/GO/KEGG)    │
│  • Generate beeswarm, bar, waterfall, heatmap, dependence     │
│  OUTPUT: disease_signature.csv, UP/DOWN gene lists            │
└──────────────────────┬───────────────────────────────────────┘
                       │
                       ▼
┌──────────────────────────────────────────────────────────────┐
│  STEP 6: Drug Repurposing  (06_drug_repurposing.py)          │
│  • Level 1: Signature reversal (CMap/L1000, Tau ≤ -90)       │
│  • Level 2: Pathway enrichment (GSEA, Subramanian PNAS 2005) │
│  • Level 3: Network pharmacology (PPI via STRING/DrugBank)    │
│  • Level 4: [Optional] Molecular docking validation           │
│  • Multi-criteria prioritization scoring                      │
│  • [Advanced] Cell-type-specific repurposing (ASGARD)         │
│  OUTPUT: repurposed_drugs.csv, drug_repurposing_summary.json  │
└──────────────────────┬───────────────────────────────────────┘
                       │
                       ▼
┌──────────────────────────────────────────────────────────────┐
│  STEP 7: Manuscript Figures  (generate_manuscript_figures.py) │
│  • Publication-quality multi-panel figures (300 DPI)           │
│  • Fig 1: Data overview & QC                                  │
│  • Fig 2: Dimensionality reduction (PCA, t-SNE)               │
│  • Fig 3: Model performance (ROC, confusion, CV)              │
│  • Fig 4: SHAP disease signature                              │
│  • Fig 5: Expression heatmap                                  │
│  • Fig 6: Drug repurposing results                            │
│  OUTPUT: figures/manuscript/*.png + *.pdf                      │
└──────────────────────────────────────────────────────────────┘
```

## Quick Start — New Project

### 1. Install Dependencies

```bash
pip install scanpy scikit-learn xgboost lightgbm shap optuna matplotlib seaborn pandas numpy scipy requests gseapy joblib
```

### 2. Project Directory Structure

```
my_project/
├── 01_load_and_qc.py
├── 02_prepare_ml_data.py
├── 03_ensemble_ml.py
├── 04_model_improvement.py
├── 05_shap_analysis.py
├── 06_drug_repurposing.py
├── generate_manuscript_figures.py
├── data/
│   ├── raw_counts.csv          # Expression matrix
│   ├── metadata.csv            # Cell annotations
│   └── ml_ready/               # Auto-generated
├── models/                     # Auto-generated
├── results/
│   └── figures/                # Auto-generated
└── manuscript/
```

### 3. Configuration (adapt per project)

Every script has a **Configuration** block at the top. The key parameters to customize:

```python
# === QC Thresholds (Step 01) ===
MIN_GENES_PER_CELL = 200       # Minimum genes detected per cell
MIN_CELLS_PER_GENE = 3         # Minimum cells expressing a gene
MAX_PCT_MITO = 20              # Max mitochondrial %. Use 5 for PBMCs, 20 for tissues

# === Label Mapping (Step 02) ===
LABEL_MAP = {
    "Disease_Group_A": 1,       # Positive class (disease)
    "Disease_Group_B": 1,       # Can merge multiple groups
    "Control_Group": 0,         # Negative class (control)
}
EXCLUDED_GROUPS = ["Confounding_Group"]  # Exclude from analysis

# === ML (Step 03) ===
TEST_SIZE = 0.20               # 80/20 train/test split
N_CV_FOLDS = 5                 # Cross-validation folds
RANDOM_STATE = 42              # Reproducibility seed
N_JOBS = -1                    # CPU cores (-1 = all)

# === SHAP (Step 05) ===
TOP_N_GENES = 100              # Disease signature size
SHAP_SAMPLE_SIZE = 2000        # Subsample for large datasets
MAX_DISPLAY = 30               # Genes shown in plots

# === Drug Repurposing (Step 06) ===
# Uses GSEApy to query L1000_Ligand_Pert databases
# Alternatively, submit UP/DOWN gene lists to clue.io manually
```

### 4. Run the Pipeline

```bash
python 01_load_and_qc.py
python 02_prepare_ml_data.py
python 03_ensemble_ml.py
python 04_model_improvement.py
python 05_shap_analysis.py
python 06_drug_repurposing.py
python generate_manuscript_figures.py
```

Each step checks for cached results and skips completed work unless `FORCE_RERUN = True`.

## Step-by-Step Instructions

### Step 1: Data Loading & Quality Control

**Goal:** Load raw scRNA-seq data, perform tiered QC, and save a processed AnnData object.

**Standard decisions:**
- **SCP format handling:** The Broad SCP metadata has a `TYPE` row as the first data row — detect and remove it
- **Matrix orientation:** Auto-detect whether genes are rows or columns by checking barcode overlap with metadata
- **MT threshold:** Use 20% for solid tissues (lung, liver), 5% for PBMCs/blood (Heumos et al., *Nat Rev Genet* 2023)
- **HVG selection:** `min_mean=0.0125, max_mean=3, min_disp=0.5` — Scanpy defaults
- **Store `adata.raw`** before HVG subsetting so all genes remain accessible for DE/SHAP

**Advanced options (recommended for publication):**
- **Ambient RNA removal:** CellBender (Fleming et al., *Nat Methods* 2023) or SoupX (Young & Behjati, *GigaScience* 2020)
- **Doublet detection:** scDblFinder (Germain et al., *F1000Research* 2021) or Scrublet (Wolock et al., *Cell Systems* 2019)
- **Normalization:** scran pooling-deconvolution (Lun et al., *Genome Biology* 2016) or SCTransform v2 Pearson residuals (Choudhary & Satija, *Genome Biology* 2022)
- **Batch correction:** Harmony (Korsunsky et al., *Nat Methods* 2019) or scVI (Lopez et al., *Nat Methods* 2018)
- **Feature selection:** Deviance-based (Townes et al., *Genome Biology* 2019) as alternative to dispersion-based HVGs

**See:** `references/qc_preprocessing.md` for detailed tiered guidance with code and tissue-specific thresholds.

### Step 2: ML Data Preparation

**Goal:** Convert the processed AnnData into ML-ready numpy arrays with binary labels.

**Key decisions:**
- **Binary labels:** Map disease groups to 1, controls to 0
- **Exclude confounders:** Remove groups that confound the disease signal (e.g., HIV-only in TB study)
- **Stratified split:** Preserve class ratios in train/test
- **Feature matrix:** Dense float32 from HVG-filtered, log-normalized counts
- **Save gene names** separately for SHAP feature naming

**See:** `references/qc_preprocessing.md` for label mapping patterns.

### Step 3: Ensemble ML Training

**Goal:** Train 4 ensemble classifiers, cross-validate, evaluate on held-out test set.

**Models:**

| Model | Type | Imbalance Handling | Key Params |
|---|---|---|---|
| Random Forest | Bagging | `class_weight='balanced'` | `n_estimators=200, min_samples_split=5` |
| XGBoost | Boosting | `scale_pos_weight=n_neg/n_pos` | `n_estimators=200, max_depth=6, lr=0.1` |
| LightGBM | Boosting | `is_unbalance=True` | `n_estimators=200, lr=0.1` |
| Stacking | Meta-learner | Inherits from base models | `final_estimator=LogisticRegression` |

**Critical: Windows threading fix for Stacking:**
```python
# Base estimators inside Stacking MUST use n_jobs=1
# to avoid nested parallelism deadlocks on Windows
StackingClassifier(
    estimators=[
        ("rf", RandomForestClassifier(n_jobs=1, ...)),
        ("xgb", XGBClassifier(n_jobs=1, ...)),
        ("lgbm", LGBMClassifier(n_jobs=1, ...)),
    ],
    n_jobs=1,  # Stacking itself also n_jobs=1
)
```

**See:** `references/ensemble_ml.md` for full model configurations and tuning strategies.

### Step 4: Model Improvement

**Goal:** Attempt hyperparameter tuning and compare against baseline performance.

**Standard:** Pre-selected grids from Bayesian optimization for speed and stability.

**Advanced options (recommended for publication):**
- **Bayesian optimization:** Optuna with TPE sampler (Akiba et al., *KDD 2019*) — more efficient than grid/random search
- **Nested cross-validation:** Inner loop for tuning, outer loop for unbiased estimation (Cawley & Talbot, *JMLR 2010*)
- **Probability calibration:** Platt scaling or isotonic regression (Niculescu-Mizil & Caruana, *ICML 2005*) — critical because tree models produce uncalibrated probabilities
- **Statistical model comparison:** DeLong test for AUC comparison (DeLong et al., *Biometrics 1988*), McNemar's test for error rate comparison
- **Extended metrics:** MCC (Chicco & Jurman, *BMC Genomics 2020*), Brier score, PR-AUC, Cohen's κ

**See:** `references/ensemble_ml.md` for full tuning strategies and statistical tests.

### Step 5: SHAP Disease Signature Extraction

**Goal:** Use SHAP to identify the genes driving disease classification and create UP/DOWN gene lists.

**Standard decisions:**
- **TreeExplainer** for exact SHAP values on tree-based models (Lundberg & Lee, *NeurIPS 2017*)
- **StackingClassifier handling:** Extract the best base estimator since TreeExplainer cannot directly explain Stacking
- **Subsample** large test sets to 2,000 cells for computational feasibility
- **Binary classification:** Extract `shap_values[1]` for the disease (positive) class
- **Signature:** Top N genes ranked by `mean(|SHAP|)`, direction = sign of `mean(SHAP)`

**Advanced options (recommended for publication):**
- **Interventional SHAP** (`feature_perturbation="interventional"`) — uses Pearl's do-calculus to avoid misattribution from correlated genes (Lundberg et al., *Nat Mach Intell* 2020; Janzing et al., *AISTATS 2020*)
- **SHAP interaction values** — captures pairwise gene interactions driving predictions (Lundberg et al., *Nat Mach Intell* 2020)
- **Bootstrap stability assessment** — 100 bootstrap resamples to identify robust signature genes (stability > 90%)
- **SHAP clustering** — cluster cells by explanation patterns to reveal disease sub-phenotypes
- **Pathway enrichment validation** — GO/KEGG/Reactome enrichment of signature genes (Subramanian et al., *PNAS 2005*)

**Disease signature output format:**

| gene | mean_abs_shap | mean_shap | direction | bootstrap_stability | rank |
|---|---|---|---|---|---|
| S100A12 | 0.0523 | +0.0523 | UP | 0.98 | 1 |
| COL3A1 | 0.0412 | -0.0412 | DOWN | 0.95 | 2 |

**See:** `references/shap_interpretation.md` for interventional SHAP, interaction values, and stability analysis.

### Step 6: Drug Repurposing

**Goal:** Multi-level drug repurposing integrating transcriptional, pathway, and network evidence.

**4-level evidence pipeline:**

1. **Level 1 — Signature Reversal (CMap/L1000):** Query via GSEApy or clue.io. Compounds with Tau ≤ -90 strongly reverse the disease state (Subramanian et al., *Cell* 2017)
2. **Level 2 — Pathway Enrichment:** Validate signature genes via GSEA/GO/KEGG enrichment (Subramanian et al., *PNAS* 2005)
3. **Level 3 — Network Pharmacology:** PPI network analysis via STRING (Szklarczyk et al., *NAR* 2023), drug-target mapping via DrugBank/DGIdb (Barabási et al., *Nat Rev Genet* 2011)
4. **Level 4 — Molecular Docking (optional):** Validate binding of top candidates to key targets (AutoDock Vina / AlphaFold structures)

**Advanced options:**
- **Single-cell-aware repurposing:** ASGARD for cell-type-specific drug identification (He et al., *Brief Bioinform* 2023)
- **Multi-criteria prioritization:** Score candidates across reversal strength, statistical significance, PPI hub targeting, DrugBank interactions, and FDA approval status
- **Consensus ranking:** Only prioritize candidates supported by ≥2 evidence levels

**See:** `references/drug_repurposing.md` for full multi-level pipeline with code.

### Step 7: Manuscript Figures

**Goal:** Generate publication-quality multi-panel figures.

**Settings:**
```python
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 10,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'axes.spines.top': False,
    'axes.spines.right': False,
})
```

**Figure catalog:**
- **Figure 1:** Data overview (class distribution, QC violins, disease groups)
- **Figure 2:** Dimensionality reduction (PCA, scree, t-SNE)
- **Figure 3:** Model performance (ROC curves, metrics bar, confusion matrices, CV)
- **Figure 4:** Disease signature (SHAP bar, direction, volcano, functional categories)
- **Figure 5:** Expression heatmap (z-scored, top 30 genes, condition-split)
- **Figure 6:** Drug repurposing (top candidates, significance volcano, thresholds, distribution)

**See:** `references/manuscript_figures.md` for detailed figure recipes.

## Best Practices

### Data & Reproducibility
- Always set `RANDOM_STATE = 42` across all scripts
- Save intermediate results (`.h5ad`, `.npz`, `.joblib`) for caching
- Each script should check for cached outputs before re-running
- Store gene names alongside feature matrices for SHAP traceability
- Use nested CV when reporting tuned model performance (Cawley & Talbot, *JMLR 2010*)

### Class Imbalance
- Use `class_weight='balanced'` (RF), `scale_pos_weight` (XGBoost), `is_unbalance=True` (LightGBM)
- Always use **stratified** splits and **stratified** CV
- Report **comprehensive** metrics: AUC-ROC, PR-AUC, MCC, F1, Brier score, Cohen's κ (Chicco & Jurman, *BMC Genomics 2020*)
- Prefer MCC over accuracy for imbalanced datasets — it considers all four confusion matrix quadrants

### QC Rigor
- Apply ambient RNA removal (CellBender) before downstream analysis (Fleming et al., *Nat Methods 2023*)
- Run doublet detection (scDblFinder/Scrublet) — expect ~0.8% doublets per 1,000 cells captured
- Use tissue-specific MT thresholds (5% for blood, 20% for lung/liver)
- Consider scran normalization for datasets with many zeros (Lun et al., *Genome Biology 2016*)

### Windows Compatibility
- Set `n_jobs=1` for base estimators inside StackingClassifier to avoid threading deadlocks
- Use `matplotlib.use('Agg')` for headless figure generation
- Use `os.chdir(os.path.dirname(os.path.abspath(__file__)))` for stable relative paths

### SHAP Rigor
- Use **interventional** TreeExplainer when genes are correlated (Lundberg et al., *Nat Mach Intell 2020*)
- Compute SHAP interaction values to identify co-regulatory gene modules
- Bootstrap the signature (100 resamples) and report stability scores
- Validate the signature biologically via GO/KEGG/Reactome pathway enrichment
- For StackingClassifier, extract the best base estimator for SHAP

### Drug Repurposing
- Focus on compounds with **Tau ≤ -90** (top 1% signature reversers, Subramanian et al., *Cell 2017*)
- Integrate ≥2 evidence levels (CMap + pathway + network) before prioritizing candidates
- Cross-reference with DrugBank/DGIdb for known drug-target interactions
- Prioritize FDA-approved or clinical-phase compounds for translational feasibility
- All computational candidates require experimental validation (in vitro → in vivo)

## Troubleshooting

### "No overlap between count matrix and metadata cell IDs"
The count matrix orientation doesn't match the metadata. The pipeline auto-detects orientation by checking row/column overlap with metadata barcodes.

### SHAP returns 3D array
For binary classifiers, SHAP may return shape `(n_samples, n_features, 2)`. Extract the positive class: `shap_values = shap_values[:, :, 1]`

### StackingClassifier + TreeExplainer error
TreeExplainer cannot handle Stacking directly. Extract a base estimator:
```python
base_model = stacking_model.named_estimators_['lgbm']
explainer = shap.TreeExplainer(base_model)
```

### Windows threading deadlock during Stacking CV
Set `n_jobs=1` for ALL base estimators AND the StackingClassifier itself.

### GSEApy API unreachable
The Enrichr/L1000 APIs can be intermittent. Implement an offline fallback or use clue.io manual queries.

## Reference Documentation

Detailed reference files for deep dives:

| File | Topic |
|---|---|
| `references/pipeline_architecture.md` | Data flow, file I/O contracts, caching strategy |
| `references/qc_preprocessing.md` | scRNA-seq QC, Scanpy patterns, label mapping |
| `references/ensemble_ml.md` | Model configs, imbalance handling, tuning, Windows fixes |
| `references/shap_interpretation.md` | TreeExplainer, signature extraction, visualization |
| `references/drug_repurposing.md` | CMap/L1000, GSEApy queries, result interpretation |
| `references/manuscript_figures.md` | Publication settings, figure recipes, color palettes |

## Example Implementation

The `examples/tb_granuloma/` directory contains the complete canonical reference implementation for tuberculosis lung granuloma drug repurposing. Consult these scripts when building a new project.

## Key References

### Core Pipeline
1. Subramanian A, et al. (2017). A Next Generation Connectivity Map. *Cell*, 171(6):1437–1452
2. Lundberg SM, Lee S-I. (2017). A Unified Approach to Interpreting Model Predictions. *NeurIPS 2017*
3. Lundberg SM, et al. (2020). From local explanations to global understanding. *Nat Mach Intell*, 2:56–67

### scRNA-seq QC & Preprocessing
4. Heumos L, et al. (2023). Best practices for single-cell analysis across modalities. *Nat Rev Genet*, 24:550–572
5. Luecken MD, Theis FJ. (2019). Current best practices in scRNA-seq analysis. *Mol Syst Biol*, 15(6):e8746
6. Fleming SJ, et al. (2023). Unsupervised removal of systematic background noise using CellBender. *Nat Methods*, 20:1323–1335
7. Lun ATL, Bach K, Marioni JC. (2016). Pooling across cells to normalize scRNA-seq data. *Genome Biol*, 17:75
8. Hafemeister C, Satija R. (2019). Normalization and variance stabilization using regularized NB regression. *Genome Biol*, 20:296
9. Korsunsky I, et al. (2019). Fast, sensitive and accurate integration with Harmony. *Nat Methods*, 16:1289–1296
10. Lopez R, et al. (2018). Deep generative modeling for single-cell transcriptomics. *Nat Methods*, 15:1053–1058
11. Germain P-L, et al. (2021). Doublet identification using scDblFinder. *F1000Research*, 10:979

### Machine Learning
12. Breiman L. (2001). Random Forests. *Machine Learning*, 45:5–32
13. Chen T, Guestrin C. (2016). XGBoost: A Scalable Tree Boosting System. *KDD 2016* (Best Paper)
14. Ke G, et al. (2017). LightGBM: A Highly Efficient Gradient Boosting Decision Tree. *NeurIPS 2017*
15. Wolpert DH. (1992). Stacked Generalization. *Neural Networks*, 5(2):241–259
16. Akiba T, et al. (2019). Optuna: A Next-generation Hyperparameter Optimization Framework. *KDD 2019*
17. Cawley GC, Talbot NLC. (2010). On Over-fitting in Model Selection. *JMLR*, 11:2079–2107
18. Chicco D, Jurman G. (2020). The advantages of the MCC over F1 and accuracy. *BMC Genomics*, 21:6
19. DeLong ER, et al. (1988). Comparing Areas Under Correlated ROC Curves. *Biometrics*, 44:837–845

### Interpretability
20. Chen H, Lundberg SM, Lee S-I. (2022). Algorithms to estimate Shapley value feature attributions. *Nat Mach Intell*, 4:658–662
21. Janzing D, Minorics L, Blöbaum P. (2020). Feature relevance quantification: A causal problem. *AISTATS 2020*

### Drug Repurposing
22. Pushpakom S, et al. (2019). Drug repurposing: progress, challenges and recommendations. *Nat Rev Drug Discov*, 18:41–58
23. Corsello SM, et al. (2020). Discovering the anticancer potential of non-oncology drugs. *Nat Cancer*, 1:235–248
24. Lamb J, et al. (2006). The Connectivity Map. *Science*, 313(5795):1929–1935
25. Barabási A-L, et al. (2011). Network medicine. *Nat Rev Genet*, 12:56–68
26. He B, et al. (2023). ASGARD for single-cell drug repurposing. *Brief Bioinform*, 24(1):bbac593

## Limitations
- Use this skill only when the task clearly matches the scope described above.
- Do not treat the output as a substitute for environment-specific validation, testing, or expert review.
- Stop and ask for clarification if required inputs, permissions, safety boundaries, or success criteria are missing.
- Drug candidates identified computationally require experimental validation before clinical consideration.
