# Manuscript Figures — Comprehensive Guide

> [!NOTE]
> Publication-quality figure recipes for scRNA-seq drug repurposing manuscripts. Follows Nature/Cell/Science formatting conventions. All figures 300 DPI, Arial/Helvetica font, clean spines.

---

## 1. Publication Settings

```python
import matplotlib
matplotlib.use('Agg')  # Headless rendering
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
import numpy as np

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.facecolor': 'white',
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'mathtext.fontset': 'dejavusans',
    'pdf.fonttype': 42,       # Editable text in PDF (not Type 3)
    'ps.fonttype': 42,
})
```

---

## 2. Journal Figure Sizing

### Column Widths (inches)

| Journal | Single Column | 1.5 Column | Double Column |
|---|---|---|---|
| **Nature** | 3.5 in (89 mm) | 5.1 in (130 mm) | 7.1 in (180 mm) |
| **Cell** | 3.35 in (85 mm) | — | 6.85 in (174 mm) |
| **Science** | 3.5 in (89 mm) | 5.5 in (140 mm) | 7.2 in (183 mm) |
| **PNAS** | 3.42 in (87 mm) | 4.5 in (114 mm) | 7.01 in (178 mm) |

### Figure Height
- Maximum: typically 9.5 in (240 mm)
- Rule of thumb: height ≤ 1.2 × width for multi-panel figures

---

## 3. Color Palette

```python
# Disease vs. Control
TB_COLOR = '#E63946'      # Red for disease (warm, attention)
CTRL_COLOR = '#457B9D'    # Blue for control (cool, recessive)

# Model colors (4 ensemble models)
COLORS_MODELS = {
    'RandomForest': '#2196F3',      # Blue
    'XGBoost': '#4CAF50',           # Green
    'LightGBM': '#FF9800',          # Orange
    'StackingEnsemble': '#9C27B0',  # Purple
}

# Metric colors (5 performance metrics)
METRIC_COLORS = ['#264653', '#2A9D8F', '#E9C46A', '#F4A261', '#E76F51']

# Functional gene categories
CAT_COLORS = {
    'Immune/Inflammatory': '#E63946',
    'Extracellular Matrix': '#457B9D',
    'Metabolism/Stress': '#F4A261',
    'Signalling/Kinase': '#2A9D8F',
    'Other': '#999999',
}

# Color-blind-safe alternative (Okabe-Ito)
COLORBLIND_SAFE = ['#E69F00', '#56B4E9', '#009E73', '#F0E442',
                    '#0072B2', '#D55E00', '#CC79A7', '#000000']
```

---

## 4. Panel Labels (Nature Style)

```python
def add_panel_label(ax, label, x=-0.15, y=1.05, fontsize=16):
    """Add Nature-style panel label (bold lowercase letter)."""
    ax.text(x, y, label, transform=ax.transAxes,
            fontsize=fontsize, fontweight='bold', va='top', ha='right')
```

---

## 5. Multi-Panel Layout with GridSpec

```python
# Standard 2×3 multi-panel figure
fig = plt.figure(figsize=(14, 10))
gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.4, wspace=0.35)

ax_a = fig.add_subplot(gs[0, 0])
ax_b = fig.add_subplot(gs[0, 1])
ax_c = fig.add_subplot(gs[0, 2])
ax_d = fig.add_subplot(gs[1, 0])
ax_e = fig.add_subplot(gs[1, 1])
ax_f = fig.add_subplot(gs[1, 2])

# Add panel labels
for ax, label in zip([ax_a, ax_b, ax_c, ax_d, ax_e, ax_f],
                      ['a', 'b', 'c', 'd', 'e', 'f']):
    add_panel_label(ax, label)

# Sub-gridspec for nested panels (e.g., 2×2 confusion matrices in one panel)
gs_inner = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs[1, 2],
                                             hspace=0.4, wspace=0.3)
```

---

## 6. Figure Catalog & Recipes

### Figure 1: Data Overview & QC

```python
# (a) Class distribution bar chart (train/test)
ax.bar(['Train Disease', 'Train Control', 'Test Disease', 'Test Control'],
       [n_train_pos, n_train_neg, n_test_pos, n_test_neg],
       color=[TB_COLOR, CTRL_COLOR, TB_COLOR, CTRL_COLOR], alpha=0.8)

# (b) Class imbalance donut chart
wedges, texts, autotexts = ax.pie(
    [n_disease, n_control], labels=['Disease', 'Control'],
    colors=[TB_COLOR, CTRL_COLOR], autopct='%1.1f%%',
    wedgeprops={'width': 0.4, 'edgecolor': 'white'})

# (c) Original disease group distribution
ax.barh(group_names, group_counts, color=group_colors)

# (d-f) QC violin plots
for metric, ax in zip(['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], axes):
    parts = ax.violinplot([adata.obs[metric].values], showmeans=True, showmedians=True)
```

### Figure 2: Dimensionality Reduction

**Critical layering:** Plot majority class FIRST (background), then minority ON TOP for visibility:

```python
# PCA scatter (PC1 vs PC2)
# Majority first (lighter, smaller)
ax.scatter(X_pca[majority_mask, 0], X_pca[majority_mask, 1],
           c=TB_COLOR, alpha=0.15, s=2, zorder=1, label='Disease')
# Minority on top (darker, larger, with edges)
ax.scatter(X_pca[minority_mask, 0], X_pca[minority_mask, 1],
           c=CTRL_COLOR, alpha=0.6, s=8, edgecolors='black',
           linewidth=0.2, zorder=2, label='Control')
```

### Figure 3: Model Performance

```python
# (a) ROC curves (all 4 models)
for name, color in COLORS_MODELS.items():
    fpr, tpr, _ = roc_curve(y_test, y_proba[name])
    auc = roc_auc_score(y_test, y_proba[name])
    ax.plot(fpr, tpr, color=color, linewidth=2,
            label=f'{name} (AUC={auc:.3f})')
ax.plot([0, 1], [0, 1], 'k--', linewidth=0.8)
ax.set_xlabel('False Positive Rate')
ax.set_ylabel('True Positive Rate')
ax.legend(fontsize=8)

# (b) Multi-metric bar chart
metrics = ['Accuracy', 'F1', 'Precision', 'Recall', 'AUC-ROC']
x = np.arange(len(metrics))
width = 0.2
for i, (name, color) in enumerate(COLORS_MODELS.items()):
    ax.bar(x + i * width, values[name], width, color=color, label=name)

# (c) Confusion matrices (2×2 nested grid)
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
for i, (name, model) in enumerate(models.items()):
    ax_sub = fig.add_subplot(gs_inner[i // 2, i % 2])
    cm = confusion_matrix(y_test, y_pred[name])
    ConfusionMatrixDisplay(cm, display_labels=['Control', 'Disease']).plot(ax=ax_sub)
    ax_sub.set_title(name, fontsize=9)
```

### Figure 4: SHAP Disease Signature

```python
# (a) Top 20 genes by mean |SHAP| colored by direction
top20 = signature.head(20).sort_values('mean_abs_shap')
colors = ['#E63946' if d == 'UP' else '#457B9D' for d in top20['direction']]
ax.barh(range(len(top20)), top20['mean_abs_shap'], color=colors)
ax.set_yticks(range(len(top20)))
ax.set_yticklabels(top20['gene'], fontsize=9)

# (b) Signature composition pie
ax.pie([n_up, n_down], labels=[f'UP ({n_up})', f'DOWN ({n_down})'],
       colors=['#E63946', '#457B9D'], autopct='%1.0f%%')

# (c) Volcano-style (mean SHAP vs mean |SHAP|)
ax.scatter(signature['mean_shap'], signature['mean_abs_shap'],
           alpha=0.3, s=10, c='gray')
# Highlight top genes
top_mask = signature['rank'] <= 20
ax.scatter(signature.loc[top_mask, 'mean_shap'],
           signature.loc[top_mask, 'mean_abs_shap'],
           s=30, c='#E63946', zorder=3)
```

### Figure 5: Expression Heatmap

```python
import seaborn as sns

# Z-score the top 30 genes
top30_genes = signature.head(30)['gene'].tolist()
gene_idx = [gene_names.index(g) for g in top30_genes]
expr = X_test[:, gene_idx]
z_scored = (expr - expr.mean(axis=0)) / (expr.std(axis=0) + 1e-8)

# Split by condition
ctrl_z = z_scored[y_test == 0]
disease_z = z_scored[y_test == 1]
combined = np.vstack([ctrl_z, disease_z])

# Heatmap
sns.heatmap(combined.T, cmap='RdBu_r', center=0, vmin=-3, vmax=3,
            yticklabels=top30_genes, xticklabels=False, ax=ax)
# Add condition separator line
ax.axvline(x=len(ctrl_z), color='black', linewidth=2)
```

### Figure 6: Drug Repurposing Results

```python
# (a) Top 15 drug candidates
top15 = drugs.head(15).sort_values('priority_score')
ax.barh(range(len(top15)), top15['priority_score'],
        color=plt.cm.RdYlGn_r(np.linspace(0.2, 0.8, len(top15))))
ax.set_yticks(range(len(top15)))
ax.set_yticklabels(top15['drug'], fontsize=9)

# (b) Significance volcano
ax.scatter(drugs['combined_score'], -np.log10(drugs['adj_pvalue']),
           alpha=0.5, c='gray', s=15)
sig_mask = drugs['adj_pvalue'] < 0.05
ax.scatter(drugs.loc[sig_mask, 'combined_score'],
           -np.log10(drugs.loc[sig_mask, 'adj_pvalue']),
           c='#E63946', s=25, zorder=3)
```

---

## 7. Statistical Annotations

```python
def add_significance(ax, x1, x2, y, p_value, height=0.02):
    """Add significance bracket between two bars."""
    # Bracket
    ax.plot([x1, x1, x2, x2], [y, y + height, y + height, y], color='black', linewidth=1)
    
    # Stars
    if p_value < 0.001:
        text = '***'
    elif p_value < 0.01:
        text = '**'
    elif p_value < 0.05:
        text = '*'
    else:
        text = 'ns'
    
    ax.text((x1 + x2) / 2, y + height, text, ha='center', va='bottom', fontsize=10)
```

---

## 8. Output Format

```python
# Always save both PNG (for review/slides) and PDF (for manuscript)
plt.savefig("figures/manuscript/figure1.png", dpi=300, facecolor='white',
            bbox_inches='tight')
plt.savefig("figures/manuscript/figure1.pdf", facecolor='white',
            bbox_inches='tight')
plt.close()
```

---

## 9. Supplementary Figures

- **Fig S1:** 3D PCA visualization (interactive with plotly if needed)
- **Fig S2:** Top feature correlation matrix (Spearman, clustered)
- **Fig S3:** Learning curves from hyperparameter tuning
- **Fig S4:** Calibration plots (reliability diagrams for each model)
- **Fig S5:** Full SHAP beeswarm plot (all genes)
- **Fig S6:** Gene-gene SHAP interaction heatmap
- **Fig S7:** Bootstrap stability distribution
- **Fig S8:** Pathway enrichment bubble plot

---

## 10. Accessibility Best Practices

- Use **color-blind-safe palettes** (Okabe-Ito, viridis, cividis)
- Add **patterns/textures** to bars when color alone distinguishes groups
- Use **different marker shapes** in scatter plots (circle, triangle, square)
- Ensure **sufficient contrast** — test with color blindness simulators
- Include **alt text** in figure captions for screen readers
