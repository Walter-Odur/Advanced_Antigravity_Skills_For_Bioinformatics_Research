# SHAP for Single-Cell Disease Signature Extraction — Comprehensive Guide

> [!NOTE]
> This guide covers the complete theory, API, visualization, and biological application of SHAP for extracting disease gene signatures from scRNA-seq ensemble ML models. Based on: [shap.readthedocs.io](https://shap.readthedocs.io) | Lundberg et al., *Nature Machine Intelligence* 2020

---

## 1. What Is SHAP?

**SHAP (SHapley Additive exPlanations)** is a game-theoretic approach to explain the output of **any** machine learning model. It connects optimal credit allocation from cooperative game theory with local explanations, providing a unified measure of feature importance.

### Why SHAP for Disease Signatures?

| Property | Benefit for scRNA-seq |
|---|---|
| **Theoretically grounded** | Shapley values from cooperative game theory with provable fairness axioms |
| **Exact for tree models** | TreeExplainer computes exact values in polynomial time for RF/XGBoost/LightGBM |
| **Local + Global** | Explains individual cell predictions AND provides global gene importance |
| **Directional** | Reveals whether a gene drives toward disease (UP) or control (DOWN) |
| **Consistent** | If a gene's true contribution increases, its SHAP value never decreases |
| **Interaction detection** | SHAP interaction values reveal gene co-regulation modules |

### Key References

1. **Shapley LS.** (1953). A Value for n-Person Games. *Contributions to the Theory of Games*, 2:307–317
2. **Lundberg SM, Lee S-I.** (2017). A Unified Approach to Interpreting Model Predictions. *NeurIPS 2017*
3. **Lundberg SM, et al.** (2020). From local explanations to global understanding with explainable AI for trees. *Nature Machine Intelligence*, 2:56–67
4. **Chen H, Lundberg SM, Lee S-I.** (2022). Algorithms to estimate Shapley value feature attributions. *Nature Machine Intelligence*, 4:658–662
5. **Janzing D, Minorics L, Blöbaum P.** (2020). Feature relevance quantification in explainable AI: A causal problem. *AISTATS 2020*
6. **Sundararajan M, Najmi A.** (2020). The many Shapley values for model explanation. *ICML 2020*

---

## 2. Mathematical Foundation

### 2.1 The Shapley Value Formula

The Shapley value for feature *i* is the **weighted average of its marginal contributions** across all possible subsets (coalitions) of features:

```
φᵢ(v) = Σ  [ |S|! × (n - |S| - 1)! / n! ] × [ v(S ∪ {i}) - v(S) ]
       S⊆N\{i}
```

**Where:**
- **N** = set of all features (genes)
- **n** = total number of features
- **S** = a subset of features excluding feature *i*
- **v(S)** = model prediction using only features in subset *S*
- **v(S ∪ {i}) - v(S)** = **marginal contribution** of gene *i* when added to coalition *S*
- The weighting factor accounts for all possible orderings of features

### 2.2 Fairness Axioms

Shapley values are the **unique** solution satisfying all four axioms simultaneously:

| Axiom | Meaning | Biological Implication |
|---|---|---|
| **Efficiency** | All SHAP values sum to `f(x) - E[f(x)]` | Total disease signal is fully distributed across genes |
| **Symmetry** | Equal-contribution features get equal values | Functionally redundant genes get equal attribution |
| **Dummy** | Zero-contribution features get zero | Non-informative genes correctly receive zero importance |
| **Linearity** | Additive across combined models | Ensemble SHAP = average of base model SHAPs |

### 2.3 Computational Complexity

Exact Shapley values require evaluating **2ⁿ** subsets (exponential in the number of features). TreeExplainer solves this in **polynomial time** for tree-based models, making it feasible for thousands of genes.

---

## 3. SHAP Explainers — Complete Reference

### 3.1 The Universal Explainer (`shap.Explainer`)

Auto-selects the best algorithm for your model:

```python
import shap

explainer = shap.Explainer(model, X_train)
shap_values = explainer(X_test)
```

### 3.2 TreeExplainer — Primary Choice for scRNA-seq ML

**Best for:** Tree-based ensemble models (XGBoost, LightGBM, Random Forest, Gradient Boosting)
**Speed:** ⚡ Very fast — computes **exact** SHAP values in polynomial time

```python
explainer = shap.TreeExplainer(model)
shap_values = explainer.shap_values(X_test)

# Modern API (returns Explanation object)
shap_values = explainer(X_test)
```

**Key Parameters:**

| Parameter | Description | Default | Recommended |
|---|---|---|---|
| `model` | Trained tree-based model | *required* | — |
| `data` | Background dataset for interventional SHAP | `None` | `X_train[:500]` |
| `model_output` | `'raw'`, `'probability'`, or `'log_loss'` | `'raw'` | `'probability'` for calibrated output |
| `feature_perturbation` | `'interventional'` or `'tree_path_dependent'` | `'tree_path_dependent'` | **`'interventional'`** for correlated genes |

### 3.3 Interventional vs. Tree-Path-Dependent SHAP

**Lundberg et al., Nature Machine Intelligence 2020**

This is the most critical decision for scRNA-seq SHAP analysis. Gene expression data is **highly correlated** — choosing the wrong mode can misattribute importance.

#### Interventional (Causal) SHAP — RECOMMENDED

Uses a background dataset to compute E[f(X) | do(Xₛ = xₛ)], following Pearl's do-calculus. **Enforces independence** between features being explained and remaining features.

```python
explainer = shap.TreeExplainer(
    model,
    data=X_train[:500],                      # Background dataset (subsample for speed)
    feature_perturbation="interventional",    # Causal framework
    model_output="probability",              # Calibrated probability output
)
shap_values = explainer.shap_values(X_explain)
```

**Why this matters for gene expression:** S100A12 and S100A8 are co-expressed in myeloid cells. Without interventional SHAP, importance may be arbitrarily split between them or misattributed to one simply because it correlates with the true driver (**Janzing et al., AISTATS 2020**).

#### Tree-Path-Dependent (Default)

Uses the tree structure itself. Faster, doesn't require background data, but can be misleading with correlated features.

```python
explainer = shap.TreeExplainer(model)  # Default: tree_path_dependent
shap_values = explainer.shap_values(X_explain)
```

#### Decision Guide

| Criterion | Interventional | Tree-Path-Dependent |
|---|---|---|
| Correlated features (genes) | ✅ Handles correctly | ⚠️ May misattribute |
| Speed | Slower (needs background) | Fast |
| Causal interpretation | Yes (do-calculus) | No (observational) |
| Background data needed | Yes | No |
| **For final disease signature** | **✅ Yes** | Quick exploration only |

### 3.4 KernelExplainer — Model-Agnostic Fallback

Use when no specialized explainer exists (e.g., neural networks, SVMs):

```python
explainer = shap.KernelExplainer(
    model.predict_proba,
    shap.sample(X_train, 100)  # Subsample for speed
)
shap_values = explainer.shap_values(X_test[:50])  # Very slow — subsample heavily
```

> [!WARNING]
> KernelExplainer is **extremely slow** for high-dimensional gene expression data. Always use TreeExplainer for tree-based models. Only use KernelExplainer as a last resort for non-tree models.

### 3.5 LinearExplainer

For logistic regression meta-learner in Stacking:

```python
explainer = shap.LinearExplainer(logistic_model, X_train)
shap_values = explainer.shap_values(X_test)
```

### Explainer Comparison Summary

| Explainer | Model Type | Exactness | Speed | For scRNA-seq? |
|---|---|---|---|---|
| **TreeExplainer** | Tree ensembles | Exact | ⚡ Very fast | **✅ Primary** |
| **LinearExplainer** | Linear models | Exact | ⚡ Very fast | Meta-learner only |
| **KernelExplainer** | Any model | Approximate | 🐢 Very slow | Last resort |
| **DeepExplainer** | Deep learning | Approximate | 🟡 Moderate | If using neural nets |

---

## 4. Handling StackingClassifier

TreeExplainer cannot directly explain a StackingClassifier because it uses a meta-learner (LogisticRegression) on top of base model predictions.

### Strategy 1: Best Base Estimator (Recommended)

```python
if best_name == "StackingEnsemble":
    # Extract best-performing base estimator
    for name, est in model.named_estimators_.items():
        if hasattr(est, "feature_importances_"):
            base_model = est
            print(f"Using base estimator: {name}")
            break
    
    explainer = shap.TreeExplainer(
        base_model,
        data=X_train[:500],
        feature_perturbation="interventional",
    )
```

### Strategy 2: Load Best Individual Model

```python
# If Stacking is best overall, use the best individual model for SHAP
test_results = pd.read_csv("results/test_results.csv", index_col=0)
non_stacking = test_results.drop("StackingEnsemble", errors="ignore")
fallback_name = non_stacking["roc_auc"].idxmax()
base_model = joblib.load(f"models/{fallback_name}.joblib")
explainer = shap.TreeExplainer(base_model, data=X_train[:500],
                                feature_perturbation="interventional")
```

### Strategy 3: Explain Each Base Model Separately

```python
# Compare SHAP across all base models for consensus ranking
all_shap = {}
for name, est in stacking_model.named_estimators_.items():
    explainer = shap.TreeExplainer(est, data=X_train[:500],
                                    feature_perturbation="interventional")
    all_shap[name] = explainer.shap_values(X_explain)

# Consensus: genes consistently ranked high across all models
consensus_importance = np.mean([np.abs(sv).mean(axis=0) for sv in all_shap.values()], axis=0)
```

---

## 5. Subsampling Strategy

### Explanation Dataset Subsampling

```python
SHAP_SAMPLE_SIZE = 2000

if X_test.shape[0] > SHAP_SAMPLE_SIZE:
    rng = np.random.RandomState(42)
    idx = rng.choice(X_test.shape[0], SHAP_SAMPLE_SIZE, replace=False)
    X_explain = X_test[idx]
    y_explain = y_test[idx]
else:
    X_explain = X_test
    y_explain = y_test
```

**Rationale:** SHAP computation scales linearly with n_samples but polynomially with model complexity. Subsampling the explanation dataset (not the model) preserves SHAP value accuracy.

### Background Dataset Subsampling (for Interventional SHAP)

```python
# Use 500 samples for background — more is diminishing returns
background = X_train[:500]
explainer = shap.TreeExplainer(model, data=background,
                                feature_perturbation="interventional")
```

---

## 6. Binary Classification SHAP Values

For binary classifiers, SHAP returns values for both classes. **Always extract class 1 (disease):**

```python
shap_values = explainer.shap_values(X_explain)

# Case 1: List of arrays [class_0_values, class_1_values]
if isinstance(shap_values, list):
    shap_values_disease = shap_values[1]

# Case 2: 3D array (n_samples, n_features, n_classes)
elif isinstance(shap_values, np.ndarray) and shap_values.ndim == 3:
    shap_values_disease = shap_values[:, :, 1]

# Case 3: 2D array (already single-class output)
else:
    shap_values_disease = shap_values

# Also extract base value (expected value)
if isinstance(explainer.expected_value, (list, np.ndarray)):
    base_value = explainer.expected_value[1]
else:
    base_value = explainer.expected_value
```

---

## 7. Disease Signature Extraction — Complete Workflow

### 7.1 Gene Importance Ranking

```python
import numpy as np
import pandas as pd

# Mean |SHAP| = global feature importance (most robust ranking)
mean_abs_shap = np.abs(shap_values_disease).mean(axis=0)

# Mean SHAP (signed) = direction of effect
mean_shap = shap_values_disease.mean(axis=0)

# Median |SHAP| (robust to outliers — consider for noisy scRNA-seq data)
median_abs_shap = np.median(np.abs(shap_values_disease), axis=0)

# Standard deviation of SHAP (variability of effect across cells)
std_shap = np.std(shap_values_disease, axis=0)

# Feature activation rate (fraction of cells where gene has non-zero SHAP)
pct_nonzero = (np.abs(shap_values_disease) > 1e-8).mean(axis=0)
```

### 7.2 Signature DataFrame

```python
signature = pd.DataFrame({
    "gene": gene_names,
    "mean_abs_shap": mean_abs_shap,
    "mean_shap": mean_shap,
    "median_abs_shap": median_abs_shap,
    "std_shap": std_shap,
    "direction": np.where(mean_shap > 0, "UP", "DOWN"),
    "pct_nonzero": pct_nonzero,
})
signature = signature.sort_values("mean_abs_shap", ascending=False)
signature["rank"] = range(1, len(signature) + 1)

# Save full ranking
signature.to_csv("results/full_gene_ranking.csv", index=False)

# Top N disease signature
TOP_N_GENES = 100
top_genes = signature.head(TOP_N_GENES)
top_genes.to_csv("results/disease_signature.csv", index=False)
```

### 7.3 UP/DOWN Gene Lists for CMap

```python
up_genes = top_genes[top_genes["direction"] == "UP"]["gene"].tolist()
down_genes = top_genes[top_genes["direction"] == "DOWN"]["gene"].tolist()

# Save as one gene per line (CMap/clue.io input format)
with open("results/signature_UP_genes.txt", "w") as f:
    f.write("\n".join(up_genes))
with open("results/signature_DOWN_genes.txt", "w") as f:
    f.write("\n".join(down_genes))

# Summary
print(f"Disease signature: {len(up_genes)} UP + {len(down_genes)} DOWN = {TOP_N_GENES} genes")
```

---

## 8. Signature Stability Assessment via Bootstrap

Assess how robust the signature is to sampling variability:

```python
n_bootstraps = 100
rng = np.random.RandomState(42)
bootstrap_rankings = []

for b in range(n_bootstraps):
    # Resample with replacement
    idx = rng.choice(shap_values_disease.shape[0], shap_values_disease.shape[0], replace=True)
    boot_importance = np.abs(shap_values_disease[idx]).mean(axis=0)
    boot_rank = np.argsort(-boot_importance)  # Descending order
    bootstrap_rankings.append(boot_rank)

# Stability: fraction of bootstraps where gene appears in top-K
K = TOP_N_GENES
stability = np.zeros(len(gene_names))
for ranking in bootstrap_rankings:
    top_k_set = set(ranking[:K])
    for gene_idx in range(len(gene_names)):
        if gene_idx in top_k_set:
            stability[gene_idx] += 1
stability /= n_bootstraps

# Add to signature
signature["bootstrap_stability"] = stability

# Report
n_stable = (signature.head(K)["bootstrap_stability"] > 0.90).sum()
print(f"Signature stability: {n_stable}/{K} genes appear in top-{K} in >90% of bootstraps")

# Flag unstable genes
unstable = signature.head(K)[signature.head(K)["bootstrap_stability"] < 0.50]
if len(unstable) > 0:
    print(f"⚠️ {len(unstable)} unstable genes (stability < 50%): {unstable['gene'].tolist()}")
```

---

## 9. SHAP Interaction Values

**Lundberg et al., Nature Machine Intelligence 2020**

Captures pairwise interaction effects between genes — how the effect of gene A changes depending on the value of gene B. This reveals **gene co-regulation modules** driving disease classification.

```python
# Compute interaction values (n_samples × n_features × n_features)
# WARNING: This is computationally expensive for large feature sets
# Subsample to ~500 cells if needed
X_interact = X_explain[:500] if X_explain.shape[0] > 500 else X_explain

shap_interaction = explainer.shap_interaction_values(X_interact)

# For binary classification, extract class 1
if isinstance(shap_interaction, list):
    shap_interaction = shap_interaction[1]

# Main effects are on the diagonal
main_effects = np.array([shap_interaction[i].diagonal() for i in range(len(X_interact))])

# Top interaction pairs (off-diagonal)
n_features = shap_interaction.shape[1]
interaction_strength = np.zeros((n_features, n_features))
for i in range(n_features):
    for j in range(n_features):
        if i != j:
            interaction_strength[i, j] = np.abs(shap_interaction[:, i, j]).mean()

# Find top interacting gene pairs
top_pairs = []
for i in range(n_features):
    for j in range(i + 1, n_features):
        top_pairs.append((gene_names[i], gene_names[j], interaction_strength[i, j]))
top_pairs.sort(key=lambda x: x[2], reverse=True)

print("Top 10 gene interactions:")
for g1, g2, strength in top_pairs[:10]:
    print(f"  {g1} × {g2}: interaction strength = {strength:.4f}")
```

### Visualizing Gene Interactions

```python
# Interaction dependence plot
shap.dependence_plot(
    (gene_names.index("S100A12"), gene_names.index("GNLY")),
    shap_interaction, X_interact,
    feature_names=gene_names, show=False
)
plt.savefig("results/figures/shap/interaction_S100A12_GNLY.png", dpi=200, bbox_inches="tight")
plt.close("all")
```

**Biological application:** Identifies gene regulatory interactions — e.g., "S100A12 SHAP value is high only when GNLY is also highly expressed" reveals co-regulatory modules in the innate immune response.

---

## 10. Pathway Enrichment of SHAP Signature

**Subramanian et al., PNAS 2005** (GSEA); **Ashburner et al., Nature Genetics 2000** (GO)

Validate the disease signature biologically using pathway enrichment:

### Over-Representation Analysis (ORA)

```python
import gseapy as gp

# Gene Ontology enrichment of top signature genes
go_results = gp.enrichr(
    gene_list=top_genes["gene"].tolist(),
    gene_sets=[
        'GO_Biological_Process_2023',
        'GO_Molecular_Function_2023',
        'GO_Cellular_Component_2023',
        'KEGG_2021_Human',
        'Reactome_2022',
        'WikiPathways_2023_Human',
        'MSigDB_Hallmark_2020',
    ],
    organism='Human',
    outdir='results/pathway_enrichment',
)

# Display top enriched pathways
top_pathways = go_results.results.sort_values('Adjusted P-value').head(20)
print(top_pathways[['Term', 'Adjusted P-value', 'Combined Score', 'Genes']].to_string())
```

### Pre-Ranked GSEA (uses full SHAP ranking, not just top-K)

```python
# Use signed SHAP values as ranking metric
rnk = signature[['gene', 'mean_shap']].copy()
rnk.columns = ['Gene', 'Score']

gsea_results = gp.prerank(
    rnk=rnk,
    gene_sets=['MSigDB_Hallmark_2020', 'KEGG_2021_Human', 'Reactome_2022'],
    outdir='results/gsea_prerank',
    min_size=15,
    max_size=500,
    permutation_num=1000,
    seed=42,
)
```

### Separate UP and DOWN Pathway Analysis

```python
# Pathways enriched in UP-regulated disease genes
up_pathways = gp.enrichr(
    gene_list=up_genes,
    gene_sets=['GO_Biological_Process_2023', 'KEGG_2021_Human'],
    organism='Human',
    outdir='results/enrichr_up_genes',
)

# Pathways enriched in DOWN-regulated disease genes
down_pathways = gp.enrichr(
    gene_list=down_genes,
    gene_sets=['GO_Biological_Process_2023', 'KEGG_2021_Human'],
    organism='Human',
    outdir='results/enrichr_down_genes',
)
```

---

## 11. SHAP Visualizations — Complete Gallery

### 11.1 Beeswarm Plot (Summary Plot)

**Scope:** Global | **Purpose:** Overview of gene importance, direction, and density

Each dot = one cell. X-axis = SHAP value. Color = gene expression (red = high, blue = low).

```python
# Modern API
shap.plots.beeswarm(shap.Explanation(
    values=shap_values_disease,
    base_values=np.full(len(X_explain), base_value),
    data=X_explain,
    feature_names=gene_names,
), max_display=30, show=False)
plt.savefig("results/figures/shap/shap_beeswarm.png", dpi=200, bbox_inches="tight")
plt.close("all")

# Legacy API (simpler)
shap.summary_plot(shap_values_disease, X_explain, feature_names=gene_names,
                  max_display=30, show=False)
plt.savefig("results/figures/shap/shap_beeswarm_legacy.png", dpi=200, bbox_inches="tight")
plt.close("all")
```

**Interpretation:**
- Features sorted top → bottom by importance
- Right-side dots → gene pushes prediction toward **disease**
- Left-side dots → gene pushes prediction toward **control**
- Red dots on right = high expression → disease (UP-regulated disease gene)
- Blue dots on right = low expression → disease (inverse relationship)

### 11.2 Bar Plot (Mean |SHAP|)

```python
shap.summary_plot(shap_values_disease, X_explain, feature_names=gene_names,
                  plot_type="bar", max_display=30, show=False)
plt.savefig("results/figures/shap/shap_bar.png", dpi=200, bbox_inches="tight")
plt.close("all")
```

### 11.3 Waterfall Plot (Single Cell Explanation)

```python
# Explain why a specific cell was classified as disease
cell_idx = 0
shap.plots.waterfall(shap.Explanation(
    values=shap_values_disease[cell_idx],
    base_values=base_value,
    data=X_explain[cell_idx],
    feature_names=gene_names,
), max_display=15, show=False)
plt.savefig("results/figures/shap/shap_waterfall_cell0.png", dpi=200, bbox_inches="tight")
plt.close("all")
```

### 11.4 Force Plot

```python
shap.initjs()

# Single cell explanation
shap.force_plot(base_value, shap_values_disease[0], X_explain[0],
                feature_names=gene_names, matplotlib=True, show=False)
plt.savefig("results/figures/shap/shap_force_cell0.png", dpi=200, bbox_inches="tight")
plt.close("all")
```

### 11.5 Dependence Plots (Top Genes)

```python
# Automatic interaction coloring for top 5 genes
top_5_genes = signature.head(5)["gene"].tolist()

for gene_name in top_5_genes:
    col_idx = gene_names.index(gene_name)
    shap.dependence_plot(col_idx, shap_values_disease, X_explain,
                         feature_names=gene_names, show=False)
    plt.savefig(f"results/figures/shap/dependence_{gene_name}.png",
                dpi=200, bbox_inches="tight")
    plt.close("all")
```

### 11.6 Heatmap Plot (Cell-Level Patterns)

```python
# Show SHAP patterns across cells — reveals disease sub-phenotypes
shap.plots.heatmap(shap.Explanation(
    values=shap_values_disease[:100],
    base_values=np.full(100, base_value),
    data=X_explain[:100],
    feature_names=gene_names,
), max_display=20, show=False)
plt.savefig("results/figures/shap/shap_heatmap.png", dpi=200, bbox_inches="tight")
plt.close("all")
```

### 11.7 Decision Plot

```python
# Shows cumulative SHAP values as genes are added
shap.decision_plot(base_value, shap_values_disease[:50],
                   feature_names=gene_names, show=False)
plt.savefig("results/figures/shap/shap_decision.png", dpi=200, bbox_inches="tight")
plt.close("all")
```

### 11.8 Signature Direction Plot (Custom)

```python
import matplotlib.pyplot as plt

top30 = signature.head(30).sort_values("mean_shap")
colors = ["#E63946" if d == "UP" else "#457B9D" for d in top30["direction"]]

fig, ax = plt.subplots(figsize=(8, 10))
ax.barh(range(len(top30)), top30["mean_shap"], color=colors, edgecolor='white', linewidth=0.5)
ax.set_yticks(range(len(top30)))
ax.set_yticklabels(top30["gene"], fontsize=9)
ax.axvline(x=0, color="black", linewidth=0.8)
ax.set_xlabel("Mean SHAP Value (disease direction →)")
ax.set_title("Disease Signature: Gene Direction")
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Legend
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color="#E63946", label="UP in disease"),
                    Patch(color="#457B9D", label="DOWN in disease")],
          loc="lower right")
plt.tight_layout()
plt.savefig("results/figures/shap/signature_direction.png", dpi=200, bbox_inches="tight")
plt.close("all")
```

### Visualization Quick Reference

| Plot | Scope | Best For | API |
|---|---|---|---|
| **Beeswarm** | Global | Gene importance + direction per cell | `shap.plots.beeswarm()` |
| **Bar** | Global | Gene ranking by mean importance | `shap.summary_plot(plot_type="bar")` |
| **Waterfall** | Local | Step-by-step single cell explanation | `shap.plots.waterfall()` |
| **Force** | Local | Compact cell prediction explanation | `shap.force_plot()` |
| **Dependence** | Global | Gene expression vs. SHAP relationship | `shap.dependence_plot()` |
| **Heatmap** | Global | Multi-cell SHAP overview | `shap.plots.heatmap()` |
| **Decision** | Both | Cumulative gene effects | `shap.decision_plot()` |
| **Direction** | Global | UP/DOWN signature composition | Custom matplotlib |

---

## 12. SHAP Clustering for Disease Sub-Phenotype Discovery

**Lundberg et al., Nature Machine Intelligence 2020**

Cluster cells by their SHAP explanation patterns (not raw expression) to discover distinct disease sub-phenotypes:

```python
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

# Cluster on SHAP values (not raw features)
# Try different numbers of clusters
best_k, best_score = 2, -1
for k in range(2, 8):
    kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
    labels = kmeans.fit_predict(shap_values_disease)
    score = silhouette_score(shap_values_disease, labels)
    if score > best_score:
        best_k, best_score = k, score

print(f"Optimal clusters: {best_k} (silhouette score: {best_score:.3f})")

# Final clustering
kmeans = KMeans(n_clusters=best_k, random_state=42, n_init=10)
shap_clusters = kmeans.fit_predict(shap_values_disease)

# Analyze per-cluster SHAP profiles
for c in range(best_k):
    mask = shap_clusters == c
    cluster_importance = np.abs(shap_values_disease[mask]).mean(axis=0)
    top_idx = np.argsort(-cluster_importance)[:5]
    print(f"\nCluster {c} ({mask.sum()} cells):")
    for i in top_idx:
        print(f"  {gene_names[i]}: mean |SHAP| = {cluster_importance[i]:.4f}")
```

**Biological insight:** Different disease sub-phenotypes may be driven by different gene modules. Cluster 0 might be inflammation-driven (S100A12, IL1B high SHAP), while Cluster 1 might be fibrosis-driven (COL3A1, FN1 high SHAP).

---

## 13. SHAP vs. Alternative Interpretability Methods

| Method | Type | Speed | Consistency | Recommended? | Reference |
|---|---|---|---|---|---|
| **SHAP (TreeExplainer)** | Exact Shapley values | Fast (polynomial) | ✅ Provably consistent | **✅ Primary** | Lundberg et al., NMI 2020 |
| Gini/Impurity importance | Heuristic (trees only) | Very fast | ❌ Biased toward high-cardinality | No | Breiman 2001 |
| Permutation importance | Model-agnostic | Slow | ⚠️ Affected by correlations | Secondary validation | Breiman 2001; Altmann et al. 2010 |
| LIME | Local surrogate | Moderate | ❌ Unstable, inconsistent | No | Ribeiro et al. 2016 |
| Integrated Gradients | Gradient-based | Fast | ⚠️ Only for differentiable models | For deep learning only | Sundararajan et al. 2017 |
| DeepLIFT | Reference-based | Fast | ⚠️ Only for neural networks | For deep learning only | Shrikumar et al. 2017 |

**For scRNA-seq ensemble ML:** Always use TreeExplainer as the primary method. Use permutation importance as a secondary validation to check that SHAP rankings are consistent.

---

## 14. Interpretation Guide

### SHAP Value Interpretation

| SHAP Value | Meaning | Gene Interpretation |
|---|---|---|
| **Positive** | Pushes prediction toward disease (class 1) | Gene contributes to disease classification |
| **Negative** | Pushes prediction toward control (class 0) | Gene contributes to control classification |
| **Large magnitude** | Strong influence on prediction | High-confidence disease gene |
| **Small magnitude** | Weak influence on prediction | Low-importance gene |

### Direction Classification

| Direction | Criterion | Biological Meaning |
|---|---|---|
| **UP** | mean_shap > 0 | Gene is upregulated in disease — higher expression → more disease-like |
| **DOWN** | mean_shap < 0 | Gene is downregulated in disease — higher expression → more control-like |

### For CMap Drug Repurposing

- **UP genes** → Find drugs that **downregulate** them (signature reversers)
- **DOWN genes** → Find drugs that **upregulate** them (signature reversers)
- A drug with **negative Tau** against the UP gene list = potential therapeutic

---

## 15. Common Pitfalls & Troubleshooting

| Pitfall | Problem | Solution |
|---|---|---|
| **Correlated genes** | SHAP splits importance between co-expressed genes | Use `feature_perturbation='interventional'` |
| **Confusing importance with causality** | High SHAP ≠ causal; reflects model reliance | Validate with pathway enrichment + literature |
| **SHAP returns 3D array** | Binary classifier outputs shape (n, p, 2) | Extract class 1: `shap_values[:, :, 1]` |
| **StackingClassifier error** | TreeExplainer can't handle meta-learner | Extract base estimator (Section 4) |
| **Interpreting overfitted model** | SHAP explains a bad model faithfully | Validate model first (AUC, calibration) |
| **Too slow** | Large test set + many genes | Subsample to 2,000 cells |
| **Negative stability** | Genes unstable across bootstraps | Report stability scores; flag genes < 50% |
| **check_additivity error** | Numerical precision issues | Set `check_additivity=False` |

### Save SHAP Values for Reuse

```python
# SHAP values are expensive to compute — always cache them
np.save("results/shap_values_disease.npy", shap_values_disease)
np.save("results/shap_base_value.npy", base_value)

# Reload later
shap_values_disease = np.load("results/shap_values_disease.npy")
base_value = np.load("results/shap_base_value.npy")
```

---

## 16. Complete End-to-End SHAP Workflow

```python
import shap
import numpy as np
import pandas as pd
import joblib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ── 1. Load model and data ──
best_model = joblib.load("models/LightGBM.joblib")
data = np.load("data/ml_ready/test_data.npz")
X_test, y_test = data["X"], data["y"]
gene_names = pd.read_csv("data/ml_ready/gene_names.csv")["gene"].tolist()

# ── 2. Subsample for computation ──
SHAP_SAMPLE_SIZE = 2000
rng = np.random.RandomState(42)
if X_test.shape[0] > SHAP_SAMPLE_SIZE:
    idx = rng.choice(X_test.shape[0], SHAP_SAMPLE_SIZE, replace=False)
    X_explain = X_test[idx]
else:
    X_explain = X_test

# ── 3. Compute interventional SHAP ──
train_data = np.load("data/ml_ready/train_data.npz")
background = train_data["X"][:500]

explainer = shap.TreeExplainer(
    best_model, data=background,
    feature_perturbation="interventional",
)
shap_values = explainer.shap_values(X_explain)

# Extract disease class
if isinstance(shap_values, list):
    shap_values_disease = shap_values[1]
elif shap_values.ndim == 3:
    shap_values_disease = shap_values[:, :, 1]
else:
    shap_values_disease = shap_values

# ── 4. Build signature ──
signature = pd.DataFrame({
    "gene": gene_names,
    "mean_abs_shap": np.abs(shap_values_disease).mean(axis=0),
    "mean_shap": shap_values_disease.mean(axis=0),
    "direction": np.where(shap_values_disease.mean(axis=0) > 0, "UP", "DOWN"),
})
signature = signature.sort_values("mean_abs_shap", ascending=False)
signature["rank"] = range(1, len(signature) + 1)

# ── 5. Bootstrap stability ──
# ... (see Section 8)

# ── 6. Save outputs ──
signature.to_csv("results/full_gene_ranking.csv", index=False)
top100 = signature.head(100)
top100.to_csv("results/disease_signature.csv", index=False)

up = top100[top100["direction"] == "UP"]["gene"].tolist()
down = top100[top100["direction"] == "DOWN"]["gene"].tolist()
with open("results/signature_UP_genes.txt", "w") as f: f.write("\n".join(up))
with open("results/signature_DOWN_genes.txt", "w") as f: f.write("\n".join(down))

# ── 7. Visualizations ──
shap.summary_plot(shap_values_disease, X_explain, feature_names=gene_names,
                  max_display=30, show=False)
plt.savefig("results/figures/shap/beeswarm.png", dpi=200, bbox_inches="tight")
plt.close("all")

# ── 8. Pathway enrichment validation ──
# ... (see Section 10)

print(f"✅ Disease signature: {len(up)} UP + {len(down)} DOWN = 100 genes")
```

---

## 17. Key References

1. Shapley LS. (1953). A Value for n-Person Games. *Contributions to the Theory of Games*, 2:307–317
2. Lundberg SM, Lee S-I. (2017). A Unified Approach to Interpreting Model Predictions. *NeurIPS 2017*
3. Lundberg SM, et al. (2020). From local explanations to global understanding. *Nat Mach Intell*, 2:56–67
4. Chen H, Lundberg SM, Lee S-I. (2022). Algorithms to estimate Shapley value feature attributions. *Nat Mach Intell*, 4:658–662
5. Janzing D, Minorics L, Blöbaum P. (2020). Feature relevance quantification: A causal problem. *AISTATS 2020*
6. Sundararajan M, Najmi A. (2020). The many Shapley values for model explanation. *ICML 2020*
7. Subramanian A, et al. (2005). Gene set enrichment analysis. *PNAS*, 102(43):15545–15550
8. Altmann A, et al. (2010). Permutation importance: a corrected feature importance measure. *Bioinformatics*, 26(10):1340–1347
