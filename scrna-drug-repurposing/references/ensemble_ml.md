# Ensemble ML Training Reference

## Literature-Grounded Best Practices

- **Breiman (2001)** — "Random Forests." *Machine Learning*, 45:5–32
- **Chen & Guestrin (2016)** — "XGBoost: A Scalable Tree Boosting System." *KDD 2016*, ACM (Best Paper Award)
- **Ke et al. (2017)** — "LightGBM: A Highly Efficient Gradient Boosting Decision Tree." *NeurIPS 2017*
- **Wolpert (1992)** — "Stacked Generalization." *Neural Networks*, 5(2):241–259
- **Lundberg et al. (2020)** — "From local explanations to global understanding with explainable AI for trees." *Nature Machine Intelligence*, 2:56–67
- **DeLong et al. (1988)** — "Comparing the Areas Under Two or More Correlated ROC Curves." *Biometrics*, 44:837–845
- **Niculescu-Mizil & Caruana (2005)** — "Predicting Good Probabilities with Supervised Learning." *ICML 2005*
- **Platt (1999)** — "Probabilistic Outputs for SVMs and Comparisons to Regularized Likelihood Methods." *Advances in Large Margin Classifiers*

---

## Model Configurations

### Random Forest (Bagging)
**Breiman, Machine Learning 2001**

Constructs an ensemble of de-correlated decision trees via bootstrap aggregation and random feature subsets. Reduces variance without increasing bias.

```python
from sklearn.ensemble import RandomForestClassifier

rf = RandomForestClassifier(
    n_estimators=200,         # 200 trees (Breiman: more trees → lower variance)
    max_depth=None,           # Grow full trees (bagging relies on full-depth trees)
    min_samples_split=5,
    min_samples_leaf=2,
    max_features='sqrt',      # √p features per split (Breiman's recommendation)
    class_weight="balanced",  # Inversely proportional to class frequency
    random_state=42,
    n_jobs=-1,
    oob_score=True,           # Out-of-bag estimate (free validation)
)
```

**Key insight:** OOB score provides an unbiased estimate of generalization error without a separate validation set (**Breiman 2001**).

### XGBoost (Gradient Boosting)
**Chen & Guestrin, KDD 2016**

Second-order gradient boosting with L1/L2 regularization. Includes built-in handling of missing values and column subsampling.

```python
from xgboost import XGBClassifier

scale_pos_weight = n_negative / n_positive

xgb = XGBClassifier(
    n_estimators=200,
    max_depth=6,              # Shallow trees for boosting (reduces overfitting)
    learning_rate=0.1,        # Shrinkage parameter η
    min_child_weight=5,       # Minimum sum of instance weight in child (regularization)
    subsample=0.8,            # Row subsampling (stochastic GB, Friedman 2002)
    colsample_bytree=0.8,     # Feature subsampling per tree
    reg_alpha=0.1,            # L1 regularization on leaf weights
    reg_lambda=1.0,           # L2 regularization on leaf weights
    scale_pos_weight=scale_pos_weight,
    eval_metric="logloss",
    random_state=42,
    n_jobs=-1,
    verbosity=0,
)
```

### LightGBM (Gradient Boosting — Leaf-wise)
**Ke et al., NeurIPS 2017**

Uses leaf-wise tree growth (vs. level-wise in XGBoost), Gradient-based One-Side Sampling (GOSS), and Exclusive Feature Bundling (EFB) for speed.

```python
from lightgbm import LGBMClassifier

lgbm = LGBMClassifier(
    n_estimators=200,
    max_depth=-1,             # No limit (leaf-wise growth is self-regularizing)
    num_leaves=31,            # Key LightGBM parameter (default, controls complexity)
    learning_rate=0.1,
    min_child_samples=20,     # Minimum data in leaf
    subsample=0.8,            # Bagging fraction
    colsample_bytree=0.8,    # Feature fraction
    reg_alpha=0.1,            # L1 regularization
    reg_lambda=1.0,           # L2 regularization
    is_unbalance=True,        # Auto-handle class imbalance
    random_state=42,
    n_jobs=-1,
    verbose=-1,
)
```

### Stacking Ensemble (Meta-learner)
**Wolpert, Neural Networks 1992**

Combines predictions from diverse base learners using a meta-learner. The key insight: base models should be maximally diverse (different inductive biases) for stacking to outperform any individual model.

```python
from sklearn.ensemble import StackingClassifier
from sklearn.linear_model import LogisticRegression

# CRITICAL: n_jobs=1 for base estimators inside Stacking
# (avoids nested parallelism deadlocks on Windows)
stacking = StackingClassifier(
    estimators=[
        ("rf", RandomForestClassifier(n_estimators=200, n_jobs=1, ...)),
        ("xgb", XGBClassifier(n_estimators=200, n_jobs=1, ...)),
        ("lgbm", LGBMClassifier(n_estimators=200, n_jobs=1, ...)),
    ],
    final_estimator=LogisticRegression(
        class_weight="balanced",
        max_iter=1000,
        random_state=42,
    ),
    cv=StratifiedKFold(n_splits=3, shuffle=True, random_state=42),
    stack_method="predict_proba",
    n_jobs=1,
)
```

---

## Class Imbalance Handling

| Model | Parameter | Mechanism | Reference |
|---|---|---|---|
| Random Forest | `class_weight='balanced'` | Weights inversely proportional to frequency | Breiman 2001 |
| XGBoost | `scale_pos_weight=n_neg/n_pos` | Scales positive gradient contribution | Chen & Guestrin 2016 |
| LightGBM | `is_unbalance=True` | Automatic weight adjustment | Ke et al. 2017 |
| Logistic Regression | `class_weight='balanced'` | Weighted cross-entropy loss | — |

**Advanced option — SMOTE variants** (use with caution):
```python
from imblearn.over_sampling import SMOTENC, ADASYN

# ADASYN focuses on hard-to-learn minority samples
# (He et al., IEEE IJCNN 2008)
adasyn = ADASYN(random_state=42)
X_resampled, y_resampled = adasyn.fit_resample(X_train, y_train)
```

**Caution:** Synthetic oversampling can introduce artifacts in gene expression data. Class weighting (built into the models) is generally preferred for scRNA-seq (**Heumos et al., Nat Rev Genet 2023**).

---

## Cross-Validation Strategy

### Standard: Stratified K-Fold
```python
from sklearn.model_selection import StratifiedKFold, cross_validate

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

scoring = {
    "accuracy": "accuracy",
    "f1": "f1",
    "precision": "precision",
    "recall": "recall",
    "roc_auc": "roc_auc",
    "average_precision": "average_precision",  # PR-AUC (better for imbalanced)
}
```

### Advanced: Nested Cross-Validation
**Cawley & Talbot, JMLR 2010** — "On Over-fitting in Model Selection"

Nested CV provides **unbiased** performance estimates when hyperparameter tuning is involved. Standard CV inflates performance because the same data is used for both tuning and evaluation.

```python
from sklearn.model_selection import cross_val_score, StratifiedKFold, GridSearchCV

# Outer loop: unbiased performance estimation
outer_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

# Inner loop: hyperparameter tuning
inner_cv = StratifiedKFold(n_splits=3, shuffle=True, random_state=42)

param_grid = {
    'n_estimators': [100, 200, 300],
    'max_depth': [6, 10, None],
    'learning_rate': [0.01, 0.05, 0.1],
}

# GridSearch inside nested CV
grid_search = GridSearchCV(
    estimator=XGBClassifier(random_state=42, n_jobs=1),
    param_grid=param_grid,
    cv=inner_cv,
    scoring='roc_auc',
    n_jobs=-1,
)

# Outer loop gives unbiased estimate
nested_scores = cross_val_score(
    grid_search, X_train, y_train,
    cv=outer_cv,
    scoring='roc_auc',
)
print(f"Nested CV AUC: {nested_scores.mean():.4f} ± {nested_scores.std():.4f}")
```

### Advanced: Bayesian Hyperparameter Optimization (Optuna)
**Akiba et al., KDD 2019** — "Optuna: A Next-generation Hyperparameter Optimization Framework"

Tree-structured Parzen Estimator (TPE) sampler for efficient Bayesian optimization. More efficient than grid/random search for high-dimensional parameter spaces.

```python
import optuna
from sklearn.model_selection import cross_val_score

def objective(trial):
    params = {
        'n_estimators': trial.suggest_int('n_estimators', 100, 500),
        'max_depth': trial.suggest_int('max_depth', 3, 15),
        'learning_rate': trial.suggest_float('learning_rate', 0.01, 0.3, log=True),
        'subsample': trial.suggest_float('subsample', 0.6, 1.0),
        'colsample_bytree': trial.suggest_float('colsample_bytree', 0.6, 1.0),
        'min_child_weight': trial.suggest_int('min_child_weight', 1, 10),
        'reg_alpha': trial.suggest_float('reg_alpha', 1e-8, 10.0, log=True),
        'reg_lambda': trial.suggest_float('reg_lambda', 1e-8, 10.0, log=True),
    }

    model = XGBClassifier(**params, random_state=42, n_jobs=1, verbosity=0)
    scores = cross_val_score(model, X_train, y_train, cv=inner_cv, scoring='roc_auc')
    return scores.mean()

study = optuna.create_study(direction='maximize', sampler=optuna.samplers.TPESampler(seed=42))
study.optimize(objective, n_trials=100, show_progress_bar=True)
```

---

## Probability Calibration
**Niculescu-Mizil & Caruana, ICML 2005**

Tree-based models produce **uncalibrated** predicted probabilities. Calibration ensures predicted probabilities match empirical frequencies.

### Calibration Assessment
```python
from sklearn.calibration import calibration_curve, CalibratedClassifierCV
from sklearn.metrics import brier_score_loss

# Brier score (lower = better calibrated)
brier = brier_score_loss(y_test, y_proba)

# Calibration curve (reliability diagram)
fraction_of_positives, mean_predicted_value = calibration_curve(
    y_test, y_proba, n_bins=10, strategy='uniform'
)
```

### Post-Hoc Calibration Methods

```python
# Platt scaling (Platt, 1999) — parametric, good for small calibration sets
calibrated_model = CalibratedClassifierCV(
    model, method='sigmoid', cv=5  # Platt scaling
)

# Isotonic regression — non-parametric, needs larger calibration set
calibrated_model = CalibratedClassifierCV(
    model, method='isotonic', cv=5
)

calibrated_model.fit(X_train, y_train)
y_proba_calibrated = calibrated_model.predict_proba(X_test)[:, 1]
```

| Method | Type | Best For | Risk |
|---|---|---|---|
| Platt scaling | Parametric (sigmoid) | Small datasets, smooth distortion | Under-flexible |
| Isotonic regression | Non-parametric | Large datasets, arbitrary distortion | Overfitting on small data |

---

## Statistical Comparison of Models

### DeLong Test for AUC Comparison
**DeLong et al., Biometrics 1988**

Tests whether two correlated AUC-ROC values are statistically significantly different.

```python
# Using scipy or dedicated implementation
from scipy import stats
import numpy as np

def delong_test(y_true, y_score_a, y_score_b):
    """
    DeLong test for comparing two correlated AUCs.
    Returns z-statistic and p-value.
    """
    from sklearn.metrics import roc_auc_score
    auc_a = roc_auc_score(y_true, y_score_a)
    auc_b = roc_auc_score(y_true, y_score_b)
    
    # Bootstrap confidence interval for AUC difference
    n_bootstraps = 10000
    rng = np.random.RandomState(42)
    auc_diffs = []
    for _ in range(n_bootstraps):
        idx = rng.randint(0, len(y_true), len(y_true))
        if len(np.unique(y_true[idx])) < 2:
            continue
        auc_a_boot = roc_auc_score(y_true[idx], y_score_a[idx])
        auc_b_boot = roc_auc_score(y_true[idx], y_score_b[idx])
        auc_diffs.append(auc_a_boot - auc_b_boot)
    
    auc_diffs = np.array(auc_diffs)
    p_value = np.mean(auc_diffs <= 0) if auc_a > auc_b else np.mean(auc_diffs >= 0)
    ci_lower, ci_upper = np.percentile(auc_diffs, [2.5, 97.5])
    
    return {
        'auc_a': auc_a, 'auc_b': auc_b,
        'auc_diff': auc_a - auc_b,
        'ci_95': (ci_lower, ci_upper),
        'p_value': 2 * min(p_value, 1 - p_value),  # Two-sided
    }
```

### McNemar's Test for Classifier Agreement
**McNemar, Psychometrika 1947**

Tests whether two classifiers make significantly different errors on the same test set.

```python
from statsmodels.stats.contingency_tables import mcnemar

# Build contingency table
# b = cases where model A is correct and B is wrong
# c = cases where model B is correct and A is wrong
correct_a = (y_pred_a == y_test)
correct_b = (y_pred_b == y_test)

b = np.sum(correct_a & ~correct_b)
c = np.sum(~correct_a & correct_b)

table = np.array([[np.sum(correct_a & correct_b), b],
                   [c, np.sum(~correct_a & ~correct_b)]])

result = mcnemar(table, exact=True)
print(f"McNemar p-value: {result.pvalue:.4f}")
```

---

## Comprehensive Metrics Suite

Beyond standard metrics, report:

```python
from sklearn.metrics import (
    matthews_corrcoef,          # MCC: best single metric for imbalanced data
    balanced_accuracy_score,    # Mean recall per class
    average_precision_score,    # PR-AUC (area under precision-recall curve)
    log_loss,                   # Negative log-likelihood (calibration-sensitive)
    brier_score_loss,           # Calibration quality
    cohen_kappa_score,          # Agreement beyond chance
)

metrics = {
    "accuracy": accuracy_score(y_test, y_pred),
    "balanced_accuracy": balanced_accuracy_score(y_test, y_pred),
    "f1": f1_score(y_test, y_pred),
    "precision": precision_score(y_test, y_pred),
    "recall": recall_score(y_test, y_pred),
    "mcc": matthews_corrcoef(y_test, y_pred),
    "roc_auc": roc_auc_score(y_test, y_proba),
    "pr_auc": average_precision_score(y_test, y_proba),
    "brier_score": brier_score_loss(y_test, y_proba),
    "log_loss": log_loss(y_test, y_proba),
    "cohen_kappa": cohen_kappa_score(y_test, y_pred),
}
```

**Key insight (Chicco & Jurman, BMC Genomics 2020):** Matthews Correlation Coefficient (MCC) is the most informative single metric for binary classification on imbalanced data — it considers all four confusion matrix quadrants.

---

## Model Serialization

```python
import joblib

joblib.dump(model, "models/LightGBM.joblib")
model = joblib.load("models/LightGBM.joblib")
```

## Windows Threading Fix

Set `n_jobs=1` for all base estimators inside StackingClassifier AND for the Stacking itself, AND for the `cross_validate()` call when evaluating Stacking. Stand-alone models can use `n_jobs=-1`.

---

## Biology-Specific ML Patterns for scRNA-seq

### Batch-Aware Train/Test Splitting

If cells come from multiple patients/samples, naive random splitting causes **data leakage** — cells from the same patient appear in both train and test sets.

```python
from sklearn.model_selection import GroupShuffleSplit

# Split by patient/sample, not by cell
groups = adata.obs['patient_id'].values
gss = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=42)
train_idx, test_idx = next(gss.split(X, y, groups))

X_train, X_test = X[train_idx], X[test_idx]
y_train, y_test = y[train_idx], y[test_idx]

# Verify: no patient appears in both sets
train_patients = set(groups[train_idx])
test_patients = set(groups[test_idx])
assert len(train_patients & test_patients) == 0, "Data leakage detected!"
```

### Cell-Type Stratification

For heterogeneous tissues, ensure all cell types are represented in train and test:

```python
# Multi-column stratification (condition + cell_type)
strat_key = adata.obs['condition'].astype(str) + '_' + adata.obs['cell_type'].astype(str)

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42, stratify=strat_key
)
```

### Pathway-Level Feature Engineering

Instead of individual gene features, aggregate genes into pathway activity scores:

```python
# Pathway scores from MSigDB Hallmark gene sets
import gseapy as gp

hallmark_sets = gp.get_library('MSigDB_Hallmark_2020')

pathway_features = np.zeros((adata.n_obs, len(hallmark_sets)))
for i, (pathway, genes) in enumerate(hallmark_sets.items()):
    genes_present = [g for g in genes if g in adata.var_names]
    if len(genes_present) >= 5:
        # Use Scanpy's score_genes for AUCell-like scoring
        sc.tl.score_genes(adata, genes_present, score_name=f'pathway_{i}')
        pathway_features[:, i] = adata.obs[f'pathway_{i}'].values

# Combine with top HVG features
X_combined = np.hstack([X_hvg, pathway_features])
```

### Handling Dropout/Zero-Inflation

scRNA-seq data has many biological and technical zeros. Strategies:

```python
# 1. Binary encoding: treat as gene detected (1) vs. not detected (0)
X_binary = (X > 0).astype(np.float32)
X_combined = np.hstack([X, X_binary])  # Doubles features but captures dropout

# 2. Zero-proportion feature: fraction of zeros per cell as a QC feature
zero_frac = (X == 0).mean(axis=1, keepdims=True)
X_with_meta = np.hstack([X, zero_frac])

# 3. Gene detection rate per cell as auxiliary feature
n_detected = (X > 0).sum(axis=1, keepdims=True)
X_with_meta = np.hstack([X, n_detected])
```

### Comparison with Traditional DE Baselines

Always compare ML-derived signatures against differential expression:

```python
# Wilcoxon rank-sum DE (standard scRNA-seq approach)
sc.tl.rank_genes_groups(adata, groupby='condition', method='wilcoxon', use_raw=True)
de_genes = sc.get.rank_genes_groups_df(adata, group='disease')
de_genes = de_genes[de_genes['pvals_adj'] < 0.05].head(TOP_N_GENES)

# Overlap with SHAP signature
shap_set = set(signature.head(TOP_N_GENES)['gene'])
de_set = set(de_genes['names'])
overlap = shap_set & de_set
jaccard = len(overlap) / len(shap_set | de_set)

print(f"SHAP vs DE overlap: {len(overlap)}/{TOP_N_GENES} genes ({jaccard:.2%} Jaccard)")
print(f"SHAP-unique genes: {shap_set - de_set}")  # Interaction-driven genes
print(f"DE-unique genes: {de_set - shap_set}")     # Marginal-only genes
```

**Key insight:** SHAP captures **non-linear interactions** and **gene-gene dependencies** that simple DE (which tests each gene independently) misses. SHAP-unique genes are often in regulatory pathways where combinatorial effects matter.

---

## Key References

1. Breiman L. (2001). Random Forests. *Machine Learning*, 45:5–32
2. Chen T, Guestrin C. (2016). XGBoost: A Scalable Tree Boosting System. *KDD 2016*
3. Ke G, et al. (2017). LightGBM: A Highly Efficient Gradient Boosting Decision Tree. *NeurIPS 2017*
4. Wolpert DH. (1992). Stacked Generalization. *Neural Networks*, 5(2):241–259
5. Akiba T, et al. (2019). Optuna: A Next-generation Hyperparameter Optimization Framework. *KDD 2019*
6. Cawley GC, Talbot NLC. (2010). On Over-fitting in Model Selection. *JMLR*, 11:2079–2107
7. Niculescu-Mizil A, Caruana R. (2005). Predicting Good Probabilities with Supervised Learning. *ICML 2005*
8. Platt JC. (1999). Probabilistic Outputs for SVMs. *Advances in Large Margin Classifiers*
9. DeLong ER, et al. (1988). Comparing Areas Under Correlated ROC Curves. *Biometrics*, 44:837–845
10. Chicco D, Jurman G. (2020). The advantages of the MCC over F1 score and accuracy. *BMC Genomics*, 21:6
