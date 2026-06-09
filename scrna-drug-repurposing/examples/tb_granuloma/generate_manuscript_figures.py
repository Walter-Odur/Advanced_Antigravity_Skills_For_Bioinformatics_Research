"""
Generate All Manuscript Figures
=================================
Creates publication-quality figures for the TB drug repurposing manuscript.
Produces both individual and composite multi-panel figures.

Usage:
    python generate_manuscript_figures.py
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import seaborn as sns
from scipy import sparse
from sklearn.decomposition import PCA
from sklearn.metrics import precision_recall_curve, average_precision_score
from sklearn.metrics import roc_curve, auc
import os
import json
import warnings
warnings.filterwarnings("ignore")

# Publication settings - Arial/Helvetica for figure labels (Nature requirement)
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
    'axes.spines.top': False,
    'axes.spines.right': False,
    'mathtext.fontset': 'dejavusans',
})

# Directories
DATA_DIR = "data"
RESULTS_DIR = "results"
FIG_DIR = os.path.join(RESULTS_DIR, "figures", "manuscript")
os.makedirs(FIG_DIR, exist_ok=True)

# Color palette
TB_COLOR = '#E63946'    # Red for TB
CTRL_COLOR = '#457B9D'  # Blue for Control
COLORS_MODELS = {
    'RandomForest': '#2196F3',
    'XGBoost': '#4CAF50',
    'LightGBM': '#FF9800',
    'StackingEnsemble': '#9C27B0'
}


def load_data():
    """Load processed AnnData and ML-ready data."""
    print("Loading AnnData...")
    adata = sc.read_h5ad(os.path.join(DATA_DIR, "processed_tb_lung.h5ad"))
    
    print("Loading ML data...")
    train_data = np.load(os.path.join(DATA_DIR, "ml_ready", "train_data.npz"), allow_pickle=True)
    test_data = np.load(os.path.join(DATA_DIR, "ml_ready", "test_data.npz"), allow_pickle=True)
    train_cells = pd.read_csv(os.path.join(DATA_DIR, "ml_ready", "train_cells.csv"))
    test_cells = pd.read_csv(os.path.join(DATA_DIR, "ml_ready", "test_cells.csv"))
    gene_names = pd.read_csv(os.path.join(DATA_DIR, "ml_ready", "gene_names.csv"))
    
    X_train = train_data['X']
    y_train = train_data['y']
    X_test = test_data['X']
    y_test = test_data['y']
    
    return adata, X_train, y_train, X_test, y_test, train_cells, test_cells, gene_names


# ================================================================
# FIGURE 1: Data Overview and Quality Control
# ================================================================
def figure1_data_overview(adata, y_train, y_test):
    """
    Figure 1: Multi-panel overview of dataset and QC.
    (a) Class distribution (train/test)
    (b) QC violin plots (n_genes, total_counts, pct_mt)
    (c) Disease status distribution across original groups
    """
    print("\n--- Generating Figure 1: Data Overview ---")
    
    fig = plt.figure(figsize=(14, 10))
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.4, wspace=0.35)
    
    # --- Panel (a): Class distribution bar chart ---
    ax_a = fig.add_subplot(gs[0, 0])
    
    train_tb = np.sum(y_train == 1)
    train_ctrl = np.sum(y_train == 0)
    test_tb = np.sum(y_test == 1)
    test_ctrl = np.sum(y_test == 0)
    
    x = np.arange(2)
    width = 0.35
    bars1 = ax_a.bar(x - width/2, [train_ctrl, train_tb], width, 
                      color=[CTRL_COLOR, TB_COLOR], alpha=0.8, label='Train', edgecolor='black', linewidth=0.5)
    bars2 = ax_a.bar(x + width/2, [test_ctrl, test_tb], width,
                      color=[CTRL_COLOR, TB_COLOR], alpha=0.5, label='Test', edgecolor='black', linewidth=0.5, hatch='///')
    
    ax_a.set_xticks(x)
    ax_a.set_xticklabels(['Control', 'TB-affected'])
    ax_a.set_ylabel('Number of cells')
    ax_a.set_title('Class distribution', fontweight='bold')
    
    # Add count labels
    for bar in bars1:
        h = bar.get_height()
        ax_a.text(bar.get_x() + bar.get_width()/2., h + 50, f'{int(h):,}', ha='center', va='bottom', fontsize=8)
    for bar in bars2:
        h = bar.get_height()
        ax_a.text(bar.get_x() + bar.get_width()/2., h + 50, f'{int(h):,}', ha='center', va='bottom', fontsize=8)
    
    ax_a.legend(loc='upper left', frameon=True, edgecolor='gray')
    ax_a.text(-0.15, 1.05, 'a', transform=ax_a.transAxes, fontsize=16, fontweight='bold', va='top')
    
    # --- Panel (b): Imbalance ratio ---
    ax_b = fig.add_subplot(gs[0, 1])
    
    total_tb = train_tb + test_tb
    total_ctrl = train_ctrl + test_ctrl
    total = total_tb + total_ctrl
    
    sizes = [total_ctrl, total_tb]
    labels_pie = [f'Control\n{total_ctrl:,} ({100*total_ctrl/total:.1f}%)', 
                  f'TB-affected\n{total_tb:,} ({100*total_tb/total:.1f}%)']
    colors_pie = [CTRL_COLOR, TB_COLOR]
    
    wedges, texts = ax_b.pie(sizes, colors=colors_pie, startangle=90,
                              wedgeprops=dict(width=0.6, edgecolor='white', linewidth=2))
    ax_b.set_title('Class imbalance', fontweight='bold')
    ax_b.legend(wedges, labels_pie, loc='center', frameon=False, fontsize=8)
    ax_b.text(-0.15, 1.05, 'b', transform=ax_b.transAxes, fontsize=16, fontweight='bold', va='top')
    
    # --- Panel (c): Original disease groups ---
    ax_c = fig.add_subplot(gs[0, 2])
    
    if 'Disease_Status' in adata.obs.columns:
        disease_counts = adata.obs['Disease_Status'].value_counts()
        group_colors = {'TB': '#E63946', 'HIVTB': '#F4A261', 'Cancer Control': '#457B9D', 'HIV': '#A8DADC'}
        bars_c = ax_c.barh(range(len(disease_counts)), disease_counts.values,
                           color=[group_colors.get(g, '#999999') for g in disease_counts.index],
                           edgecolor='black', linewidth=0.5)
        ax_c.set_yticks(range(len(disease_counts)))
        ax_c.set_yticklabels(disease_counts.index)
        ax_c.set_xlabel('Number of cells')
        ax_c.set_title('Disease status groups', fontweight='bold')
        for i, (v, idx) in enumerate(zip(disease_counts.values, disease_counts.index)):
            ax_c.text(v + 100, i, f'{v:,}', va='center', fontsize=8)
    ax_c.text(-0.15, 1.05, 'c', transform=ax_c.transAxes, fontsize=16, fontweight='bold', va='top')
    
    # --- Panel (d-f): QC violin plots ---
    # Filter to ML-used cells (TB + HIVTB + Cancer Control only)
    mask = adata.obs['Disease_Status'].isin(['TB', 'HIVTB', 'Cancer Control'])
    adata_ml = adata[mask].copy()
    adata_ml.obs['condition'] = adata_ml.obs['Disease_Status'].map(
        {'TB': 'TB-affected', 'HIVTB': 'TB-affected', 'Cancer Control': 'Control'}
    )
    
    qc_metrics = [('n_genes_by_counts', 'Genes detected per cell', 'd'),
                  ('total_counts', 'Total UMI counts per cell', 'e'),
                  ('pct_counts_mt', 'Mitochondrial fraction (%)', 'f')]
    
    for i, (metric, ylabel, panel_label) in enumerate(qc_metrics):
        ax = fig.add_subplot(gs[1, i])
        if metric in adata_ml.obs.columns:
            data_ctrl = adata_ml.obs.loc[adata_ml.obs['condition'] == 'Control', metric].values
            data_tb = adata_ml.obs.loc[adata_ml.obs['condition'] == 'TB-affected', metric].values
            
            parts = ax.violinplot([data_ctrl, data_tb], positions=[0, 1], showmeans=True, showmedians=True)
            for j, pc in enumerate(parts['bodies']):
                pc.set_facecolor([CTRL_COLOR, TB_COLOR][j])
                pc.set_alpha(0.7)
            parts['cmeans'].set_color('black')
            parts['cmedians'].set_color('gray')
            parts['cmedians'].set_linestyle('--')
            
            ax.set_xticks([0, 1])
            ax.set_xticklabels(['Control', 'TB-affected'])
            ax.set_ylabel(ylabel)
            ax.set_title(ylabel, fontweight='bold')
        ax.text(-0.15, 1.05, panel_label, transform=ax.transAxes, fontsize=16, fontweight='bold', va='top')
    
    plt.savefig(os.path.join(FIG_DIR, "figure1_data_overview.png"), dpi=300, facecolor='white')
    plt.savefig(os.path.join(FIG_DIR, "figure1_data_overview.pdf"), facecolor='white')
    plt.close()
    print("  Saved Figure 1")


# ================================================================
# FIGURE 2: Dimensionality Reduction and Clustering
# ================================================================
def figure2_dimensionality_reduction(adata, X_train, y_train, X_test, y_test):
    """
    Figure 2: Dimensionality reduction visualizations.
    (a) PCA colored by condition (Control plotted ON TOP for visibility)
    (b) PCA variance explained (scree plot)
    (c) t-SNE colored by condition (Control plotted ON TOP for visibility)
    """
    print("\n--- Generating Figure 2: Dimensionality Reduction ---")
    
    fig = plt.figure(figsize=(16, 5))
    gs = gridspec.GridSpec(1, 3, figure=fig, wspace=0.35)
    
    # --- PCA on ML data ---
    print("  Computing PCA on ML features...")
    pca = PCA(n_components=50, random_state=42)
    X_pca_train = pca.fit_transform(X_train)
    
    # Panel (a): PCA scatter colored by condition
    # CRITICAL FIX: Plot TB first (background), then Control on top for visibility
    ax_a = fig.add_subplot(gs[0, 0])
    ctrl_mask = y_train == 0
    tb_mask = y_train == 1
    
    # TB first (background, lighter alpha)
    ax_a.scatter(X_pca_train[tb_mask, 0], X_pca_train[tb_mask, 1], 
                 c=TB_COLOR, alpha=0.15, s=2, label='TB-affected', rasterized=True, zorder=1)
    # Control on top (stronger alpha, larger points, dark edge)
    ax_a.scatter(X_pca_train[ctrl_mask, 0], X_pca_train[ctrl_mask, 1], 
                 c=CTRL_COLOR, alpha=0.6, s=8, label='Control', rasterized=True, 
                 edgecolors='black', linewidth=0.2, zorder=2)
    
    ax_a.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
    ax_a.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
    ax_a.set_title('PCA of ML features', fontweight='bold')
    ax_a.legend(markerscale=3, frameon=True, edgecolor='gray')
    ax_a.text(-0.12, 1.05, 'a', transform=ax_a.transAxes, fontsize=16, fontweight='bold', va='top')
    
    # Panel (b): Variance explained (scree plot)
    ax_b = fig.add_subplot(gs[0, 1])
    cum_var = np.cumsum(pca.explained_variance_ratio_) * 100
    ax_b.bar(range(1, 21), pca.explained_variance_ratio_[:20]*100, color='#457B9D', alpha=0.7, label='Individual')
    ax_b2 = ax_b.twinx()
    ax_b2.plot(range(1, 21), cum_var[:20], 'o-', color=TB_COLOR, markersize=4, label='Cumulative')
    ax_b2.set_ylabel('Cumulative variance (%)')
    ax_b2.spines['top'].set_visible(False)
    ax_b.set_xlabel('Principal Component')
    ax_b.set_ylabel('Variance explained (%)')
    ax_b.set_title('PCA variance explained', fontweight='bold')
    
    # Combined legend
    lines1, labels1 = ax_b.get_legend_handles_labels()
    lines2, labels2 = ax_b2.get_legend_handles_labels()
    ax_b.legend(lines1 + lines2, labels1 + labels2, loc='center right', frameon=True, edgecolor='gray')
    ax_b.text(-0.12, 1.05, 'b', transform=ax_b.transAxes, fontsize=16, fontweight='bold', va='top')
    
    # Panel (c): t-SNE from PCA features
    ax_c = fig.add_subplot(gs[0, 2])
    from sklearn.manifold import TSNE
    print("  Computing t-SNE on top 20 PCs (5,000-cell subsample)...")
    subsample = min(5000, X_pca_train.shape[0])
    idx = np.random.RandomState(42).choice(X_pca_train.shape[0], subsample, replace=False)
    tsne = TSNE(n_components=2, random_state=42, perplexity=30)
    X_tsne = tsne.fit_transform(X_pca_train[idx, :20])
    y_sub = y_train[idx]
    
    # CRITICAL FIX: Plot TB first (background), then Control on top
    ax_c.scatter(X_tsne[y_sub == 1, 0], X_tsne[y_sub == 1, 1],
                 c=TB_COLOR, alpha=0.15, s=2, label='TB-affected', rasterized=True, zorder=1)
    ax_c.scatter(X_tsne[y_sub == 0, 0], X_tsne[y_sub == 0, 1],
                 c=CTRL_COLOR, alpha=0.6, s=8, label='Control', rasterized=True,
                 edgecolors='black', linewidth=0.2, zorder=2)
    
    ax_c.set_xlabel('t-SNE 1')
    ax_c.set_ylabel('t-SNE 2')
    ax_c.set_title('t-SNE of ML features', fontweight='bold')
    ax_c.legend(markerscale=3, frameon=True, edgecolor='gray')
    ax_c.text(-0.12, 1.05, 'c', transform=ax_c.transAxes, fontsize=16, fontweight='bold', va='top')
    
    plt.savefig(os.path.join(FIG_DIR, "figure2_dimensionality_reduction.png"), dpi=300, facecolor='white')
    plt.savefig(os.path.join(FIG_DIR, "figure2_dimensionality_reduction.pdf"), facecolor='white')
    plt.close()
    print("  Saved Figure 2")
    
    return pca


# ================================================================
# FIGURE 3: Model Performance Comparison (with ROC curves)
# ================================================================
def figure3_model_performance(X_train, y_train, X_test, y_test):
    """
    Figure 3: Ensemble ML model performance.
    (a) ROC curves (embedded from pipeline output)
    (b) Model comparison bar chart
    (c) Confusion matrices
    (d) Cross-validation box plots
    """
    print("\n--- Generating Figure 3: Model Performance ---")
    
    fig = plt.figure(figsize=(14, 12))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.3)
    
    # Load data
    test_results = pd.read_csv(os.path.join(RESULTS_DIR, "test_results.csv"), index_col=0)
    cv_results = pd.read_csv(os.path.join(RESULTS_DIR, "cv_results.csv"), index_col=0)
    
    models = ['RandomForest', 'XGBoost', 'LightGBM', 'StackingEnsemble']
    model_labels = ['Random Forest', 'XGBoost', 'LightGBM', 'Stacking\nEnsemble']
    
    # Panel (a): ROC Curves -- load from existing pipeline output
    ax_a = fig.add_subplot(gs[0, 0])
    
    roc_img_path = os.path.join(RESULTS_DIR, "figures", "ml", "roc_curves.png")
    if os.path.exists(roc_img_path):
        roc_img = plt.imread(roc_img_path)
        ax_a.imshow(roc_img, aspect='auto')
        ax_a.set_axis_off()
    else:
        ax_a.text(0.5, 0.5, 'ROC curves\n(not available)', transform=ax_a.transAxes,
                  ha='center', va='center', fontsize=14, color='gray')
        ax_a.set_axis_off()
    ax_a.text(-0.05, 1.05, 'a', transform=ax_a.transAxes, fontsize=16, fontweight='bold', va='top')
    
    # Panel (b): Multi-metric bar chart 
    ax_b = fig.add_subplot(gs[0, 1])
    metrics = ['accuracy', 'f1', 'precision', 'recall', 'roc_auc']
    metric_labels_list = ['Accuracy', 'F1', 'Precision', 'Recall', 'AUC-ROC']
    metric_colors = ['#264653', '#2A9D8F', '#E9C46A', '#F4A261', '#E76F51']
    
    x = np.arange(len(models))
    bar_width = 0.15
    
    for i, (metric, mlabel, color) in enumerate(zip(metrics, metric_labels_list, metric_colors)):
        values = [test_results.loc[m, metric] for m in models]
        offset = (i - 2) * bar_width
        bars = ax_b.bar(x + offset, values, bar_width, label=mlabel, color=color, 
                        edgecolor='white', linewidth=0.5)
    
    ax_b.set_xticks(x)
    ax_b.set_xticklabels(['RF', 'XGB', 'LGBM', 'Stack'], fontsize=9)
    ax_b.set_ylabel('Score')
    ax_b.set_ylim(0.85, 1.01)
    ax_b.set_title('Test set performance metrics', fontweight='bold')
    ax_b.legend(loc='lower left', frameon=True, edgecolor='gray', fontsize=7, ncol=2)
    ax_b.text(-0.12, 1.05, 'b', transform=ax_b.transAxes, fontsize=16, fontweight='bold', va='top')
    
    # Panel (c): Confusion matrix heatmaps (2x2 grid)
    ax_c = fig.add_subplot(gs[1, 0])
    # Manually create from known results
    cm_data = {
        'Random Forest': np.array([[207, 227], [31, 3344]]),
        'XGBoost': np.array([[362, 72], [246, 3129]]),
        'LightGBM': np.array([[338, 96], [115, 3260]]),
        'Stacking': np.array([[392, 42], [288, 3087]])
    }
    
    inner_gs = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs[1, 0], hspace=0.4, wspace=0.3)
    for idx, (model, cm) in enumerate(cm_data.items()):
        ax_cm = fig.add_subplot(inner_gs[idx // 2, idx % 2])
        sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', cbar=False,
                    xticklabels=['Ctrl', 'TB'], yticklabels=['Ctrl', 'TB'],
                    ax=ax_cm, annot_kws={'size': 9})
        ax_cm.set_title(model, fontsize=9, fontweight='bold')
        if idx % 2 == 0:
            ax_cm.set_ylabel('True', fontsize=8)
        if idx >= 2:
            ax_cm.set_xlabel('Predicted', fontsize=8)
    
    ax_c.set_visible(False)
    ax_c_title = fig.text(0.07, 0.48, 'c', fontsize=16, fontweight='bold', va='top')
    
    # Panel (d): CV AUC-ROC comparison with error bars
    ax_d = fig.add_subplot(gs[1, 1])
    cv_means = [cv_results.loc[m, 'roc_auc_mean'] for m in models]
    cv_stds = [cv_results.loc[m, 'roc_auc_std'] for m in models]
    
    model_colors = [COLORS_MODELS[m] for m in models]
    bars_d = ax_d.barh(range(len(models)), cv_means, xerr=cv_stds,
                        color=model_colors, alpha=0.8, edgecolor='black', linewidth=0.5,
                        capsize=5, error_kw={'linewidth': 1.5})
    ax_d.set_yticks(range(len(models)))
    ax_d.set_yticklabels(['Random Forest', 'XGBoost', 'LightGBM', 'Stacking\nEnsemble'], fontsize=9)
    ax_d.set_xlabel('AUC-ROC (5-fold CV)')
    ax_d.set_xlim(0.9, 0.98)
    ax_d.set_title('Cross-validation AUC-ROC', fontweight='bold')
    for i, (mean, std) in enumerate(zip(cv_means, cv_stds)):
        ax_d.text(mean + std + 0.002, i, f'{mean:.3f} \u00b1 {std:.3f}', va='center', fontsize=8)
    ax_d.text(-0.12, 1.05, 'd', transform=ax_d.transAxes, fontsize=16, fontweight='bold', va='top')
    
    plt.savefig(os.path.join(FIG_DIR, "figure3_model_performance.png"), dpi=300, facecolor='white')
    plt.savefig(os.path.join(FIG_DIR, "figure3_model_performance.pdf"), facecolor='white')
    plt.close()
    print("  Saved Figure 3")


# ================================================================
# FIGURE 4: SHAP Feature Importance and Disease Signature
# ================================================================
def figure4_shap_signature(gene_names):
    """
    Figure 4: SHAP-derived disease signature.
    (a) SHAP beeswarm (top 20 genes)
    (b) Signature direction bar chart
    (c) Volcano-style plot (mean |SHAP| vs direction)
    (d) Gene functional categories (computed, not hardcoded)
    """
    print("\n--- Generating Figure 4: Disease Signature ---")
    
    signature = pd.read_csv(os.path.join(RESULTS_DIR, "disease_signature.csv"))
    
    fig = plt.figure(figsize=(14, 12))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.35)
    
    # Panel (a): Top 20 genes by mean |SHAP| (horizontal bar)
    ax_a = fig.add_subplot(gs[0, 0])
    top20 = signature.head(20).copy()
    top20 = top20.iloc[::-1]  # Reverse for horizontal bar
    colors_bar = [TB_COLOR if d == 'UP' else CTRL_COLOR for d in top20['direction']]
    ax_a.barh(range(len(top20)), top20['mean_abs_shap'], color=colors_bar, 
              edgecolor='black', linewidth=0.3, alpha=0.8)
    ax_a.set_yticks(range(len(top20)))
    ax_a.set_yticklabels(top20['gene'], fontsize=8, style='italic')
    ax_a.set_xlabel('Mean |SHAP value|')
    ax_a.set_title('Top 20 genes by SHAP importance', fontweight='bold')
    legend_elements = [Patch(facecolor=TB_COLOR, label='Upregulated in TB'),
                       Patch(facecolor=CTRL_COLOR, label='Downregulated in TB')]
    ax_a.legend(handles=legend_elements, loc='lower right', frameon=True, edgecolor='gray', fontsize=8)
    ax_a.text(-0.15, 1.05, 'a', transform=ax_a.transAxes, fontsize=16, fontweight='bold', va='top')
    
    # Panel (b): Signature direction - up vs down counts
    ax_b = fig.add_subplot(gs[0, 1])
    n_up = signature[signature['direction'] == 'UP'].shape[0]
    n_down = signature[signature['direction'] == 'DOWN'].shape[0]
    
    bars_dir = ax_b.bar(['Upregulated\nin TB', 'Downregulated\nin TB'], [n_up, n_down],
                         color=[TB_COLOR, CTRL_COLOR], edgecolor='black', linewidth=0.5, alpha=0.8, width=0.5)
    ax_b.set_ylabel('Number of genes')
    ax_b.set_title('Disease signature composition', fontweight='bold')
    for bar, val in zip(bars_dir, [n_up, n_down]):
        ax_b.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.5,
                  str(val), ha='center', va='bottom', fontweight='bold', fontsize=12)
    ax_b.text(-0.12, 1.05, 'b', transform=ax_b.transAxes, fontsize=16, fontweight='bold', va='top')
    
    # Panel (c): Volcano-style plot (mean_shap vs mean_abs_shap)
    ax_c = fig.add_subplot(gs[1, 0])
    sig = signature.copy()
    colors_vol = [TB_COLOR if d == 'UP' else CTRL_COLOR for d in sig['direction']]
    ax_c.scatter(sig['mean_shap'], sig['mean_abs_shap'], c=colors_vol, alpha=0.6, s=30, edgecolor='black', linewidth=0.3)
    
    # Label top 10 with adjusted text to avoid overlaps
    from matplotlib import patheffects
    for _, row in sig.head(10).iterrows():
        ax_c.annotate(row['gene'], (row['mean_shap'], row['mean_abs_shap']),
                       fontsize=7, style='italic', ha='center', va='bottom',
                       xytext=(0, 5), textcoords='offset points',
                       path_effects=[patheffects.withStroke(linewidth=2, foreground='white')])
    
    ax_c.axvline(x=0, color='gray', linestyle='--', linewidth=0.5)
    ax_c.set_xlabel('Mean SHAP value (directional)')
    ax_c.set_ylabel('Mean |SHAP value| (importance)')
    ax_c.set_title('Feature importance vs. direction', fontweight='bold')
    ax_c.text(-0.12, 1.05, 'c', transform=ax_c.transAxes, fontsize=16, fontweight='bold', va='top')
    
    # Panel (d): Gene functional categories - PROPERLY COMPUTED from top 30
    ax_d = fig.add_subplot(gs[1, 1])
    
    # Categorize top 30 genes by biological function
    categories = {
        'Immune/Inflammatory': ['S100A12', 'CCL18', 'OSM', 'IL1RN', 'GNLY', 'IFITM2', 'MARCO', 'CCL4', 'FCGR3B', 'PLEK', 'HLA-DRB5'],
        'Extracellular Matrix': ['COL3A1', 'COL1A2', 'COL14A1', 'FBLN1', 'DCN', 'SPARC', 'FN1'],
        'Metabolism/Stress': ['MTRNR2L12', 'HMOX1', 'LIPA', 'LPL', 'PLTP', 'APOC1', 'HSPA1B', 'HBB', 'FOLR3'],
        'Signalling/Kinase': ['TAOK1', 'RGS2', 'CTSL'],
    }
    
    top30_genes = set(signature.head(30)['gene'].values)
    cat_counts = {}
    for cat, genes in categories.items():
        count = len(set(genes) & top30_genes)
        if count > 0:
            cat_counts[cat] = count
    
    # Check for uncategorized
    all_categorized = set()
    for genes in categories.values():
        all_categorized.update(genes)
    uncategorized = top30_genes - all_categorized
    if uncategorized:
        cat_counts['Other'] = len(uncategorized)
    
    cat_colors = {'Immune/Inflammatory': '#E63946', 'Extracellular Matrix': '#457B9D', 
                  'Metabolism/Stress': '#F4A261', 'Signalling/Kinase': '#2A9D8F', 'Other': '#999999'}
    
    # Use horizontal bar chart instead of pie chart for more professional look
    cat_names = list(cat_counts.keys())
    cat_vals = list(cat_counts.values())
    total_genes = sum(cat_vals)
    y_pos = range(len(cat_names))
    
    bars_cat = ax_d.barh(y_pos, cat_vals, 
                          color=[cat_colors.get(k, '#999') for k in cat_names],
                          edgecolor='black', linewidth=0.3, alpha=0.85)
    ax_d.set_yticks(y_pos)
    ax_d.set_yticklabels(cat_names, fontsize=9)
    ax_d.set_xlabel('Number of genes')
    ax_d.set_title('Functional categories (top 30 genes)', fontweight='bold')
    
    for i, (bar, val) in enumerate(zip(bars_cat, cat_vals)):
        pct = 100 * val / total_genes
        ax_d.text(bar.get_width() + 0.2, i, f'{val} ({pct:.0f}%)', va='center', fontsize=8)
    
    ax_d.text(-0.12, 1.05, 'd', transform=ax_d.transAxes, fontsize=16, fontweight='bold', va='top')
    
    plt.savefig(os.path.join(FIG_DIR, "figure4_disease_signature.png"), dpi=300, facecolor='white')
    plt.savefig(os.path.join(FIG_DIR, "figure4_disease_signature.pdf"), facecolor='white')
    plt.close()
    print("  Saved Figure 4")


# ================================================================
# FIGURE 6: Drug Repurposing Results
# ================================================================
def figure6_drug_repurposing():
    """
    Figure 6: Drug repurposing analysis.
    (a) Top drug candidates by enrichment score (CLEANED names)
    (b) Significance vs effect size scatter
    (c) Drug source/class distribution
    (d) Score distribution histogram
    """
    print("\n--- Generating Figure 6: Drug Repurposing ---")
    
    drugs = pd.read_csv(os.path.join(RESULTS_DIR, "repurposed_drugs.csv"))
    
    fig = plt.figure(figsize=(14, 10))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.4, wspace=0.35)
    
    # Panel (a): Top 15 drug candidates
    ax_a = fig.add_subplot(gs[0, 0])
    
    # Get top candidates
    if 'combined_score' in drugs.columns:
        score_col = 'combined_score'
    elif 'Combined Score' in drugs.columns:
        score_col = 'Combined Score'
    else:
        score_col = drugs.columns[1]
    
    if 'Adjusted P-value' in drugs.columns:
        pval_col = 'Adjusted P-value'
    elif 'adjusted_p_value' in drugs.columns:
        pval_col = 'adjusted_p_value'
    else:
        pval_col = None
    
    if 'Term' in drugs.columns:
        name_col = 'Term'
    elif 'term' in drugs.columns:
        name_col = 'term'
    else:
        name_col = drugs.columns[0]
    
    top_drugs = drugs.nlargest(15, score_col).copy()
    top_drugs = top_drugs.iloc[::-1]  # Reverse for horizontal bar
    
    # FIXED: Clean drug names properly - extract compound name, remove cell lines/IDs
    def clean_drug_name(name):
        """Extract just the drug compound name, removing cell line info and database IDs."""
        import re
        # Remove common suffixes like 'HL60 UP', 'HL60 DOWN', 'PC3 UP', etc.
        name = re.sub(r'\s+(HL60|PC3|MCF7|A375|A549|HCC515|HEPG2|HT29|VCAP)\s+(UP|DOWN)$', '', name, flags=re.IGNORECASE)
        # Remove database IDs like 'DB00818', 'TTD 00008407', 'CTD 00006230'
        name = re.sub(r'\s+(DB\d+|TTD\s*\d+|CTD\s*\d+)\s*', ' ', name, flags=re.IGNORECASE)
        # Remove GEO sample references like 'human GSE30903 sample 3220'
        name = re.sub(r'\s+(human|rat|mouse)\s+GSE\d+\s+sample\s+\d+', '', name, flags=re.IGNORECASE)
        # Remove trailing numeric IDs like '-2429', '-2675', ' 7305', ' 3902'
        name = re.sub(r'[\s-]+\d{3,}$', '', name)
        # Clean up extra whitespace
        name = ' '.join(name.split()).strip()
        # Capitalize first letter
        if name:
            name = name[0].upper() + name[1:]
        return name
    
    display_names = [clean_drug_name(n) for n in top_drugs[name_col]]
    
    colors_drugs = plt.cm.RdYlBu_r(np.linspace(0.2, 0.8, len(top_drugs)))
    ax_a.barh(range(len(top_drugs)), top_drugs[score_col], color=colors_drugs,
              edgecolor='black', linewidth=0.3)
    ax_a.set_yticks(range(len(top_drugs)))
    ax_a.set_yticklabels(display_names, fontsize=8)
    ax_a.set_xlabel('Combined Enrichment Score')
    ax_a.set_title('Top 15 drug candidates', fontweight='bold')
    ax_a.text(-0.2, 1.05, 'a', transform=ax_a.transAxes, fontsize=16, fontweight='bold', va='top')
    
    # Panel (b): Significance volcano plot
    ax_b = fig.add_subplot(gs[0, 1])
    if pval_col:
        neg_log_p = -np.log10(drugs[pval_col].clip(1e-15))
        scatter = ax_b.scatter(drugs[score_col], neg_log_p, 
                               c=drugs[score_col], cmap='RdYlBu_r',
                               alpha=0.6, s=30, edgecolor='black', linewidth=0.2)
        ax_b.axhline(y=-np.log10(0.05), color='red', linestyle='--', linewidth=1, label='p = 0.05')
        ax_b.set_xlabel('Combined Enrichment Score')
        ax_b.set_ylabel('$-\\log_{10}$(Adjusted $P$-value)')
        ax_b.set_title('Statistical significance vs. effect size', fontweight='bold')
        ax_b.legend(frameon=True, edgecolor='gray')
        plt.colorbar(scatter, ax=ax_b, label='Score', shrink=0.8)
    ax_b.text(-0.12, 1.05, 'b', transform=ax_b.transAxes, fontsize=16, fontweight='bold', va='top')
    
    # Panel (c): Number of significant drugs at different thresholds
    ax_c = fig.add_subplot(gs[1, 0])
    if pval_col:
        thresholds = [0.05, 0.01, 0.001, 0.0001]
        threshold_labels = ['$P$ < 0.05', '$P$ < 0.01', '$P$ < 0.001', '$P$ < 0.0001']
        counts_thresh = [drugs[drugs[pval_col] < t].shape[0] for t in thresholds]
        
        bars_c = ax_c.bar(threshold_labels, counts_thresh, 
                           color=['#A8DADC', '#457B9D', '#1D3557', '#0B1F3A'],
                           edgecolor='black', linewidth=0.5)
        ax_c.set_ylabel('Number of drug candidates')
        ax_c.set_title('Drug candidates at significance thresholds', fontweight='bold')
        for bar, val in zip(bars_c, counts_thresh):
            ax_c.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 2,
                      str(val), ha='center', va='bottom', fontweight='bold')
    ax_c.text(-0.12, 1.05, 'c', transform=ax_c.transAxes, fontsize=16, fontweight='bold', va='top')
    
    # Panel (d): Score distribution histogram
    ax_d = fig.add_subplot(gs[1, 1])
    ax_d.hist(drugs[score_col], bins=50, color='#457B9D', alpha=0.7, edgecolor='black', linewidth=0.3)
    if pval_col:
        sig_scores = drugs.loc[drugs[pval_col] < 0.05, score_col]
        ax_d.hist(sig_scores, bins=50, color=TB_COLOR, alpha=0.7, edgecolor='black', linewidth=0.3, label='$P$ < 0.05')
    ax_d.set_xlabel('Combined Enrichment Score')
    ax_d.set_ylabel('Number of compounds')
    ax_d.set_title('Distribution of enrichment scores', fontweight='bold')
    ax_d.legend(frameon=True, edgecolor='gray')
    ax_d.text(-0.12, 1.05, 'd', transform=ax_d.transAxes, fontsize=16, fontweight='bold', va='top')
    
    plt.savefig(os.path.join(FIG_DIR, "figure6_drug_repurposing.png"), dpi=300, facecolor='white')
    plt.savefig(os.path.join(FIG_DIR, "figure6_drug_repurposing.pdf"), facecolor='white')
    plt.close()
    print("  Saved Figure 6")


# ================================================================
# FIGURE 5: Expression Heatmap of Top Signature Genes
# ================================================================
def figure5_expression_heatmap(adata, gene_names):
    """
    Figure 5: Expression heatmap of top SHAP genes across conditions.
    """
    print("\n--- Generating Figure 5: Expression Heatmap ---")
    
    signature = pd.read_csv(os.path.join(RESULTS_DIR, "disease_signature.csv"))
    top_genes = signature.head(30)['gene'].values
    
    # Filter adata to ML cells
    mask = adata.obs['Disease_Status'].isin(['TB', 'HIVTB', 'Cancer Control'])
    adata_ml = adata[mask].copy()
    adata_ml.obs['condition'] = adata_ml.obs['Disease_Status'].map(
        {'TB': 'TB-affected', 'HIVTB': 'TB-affected', 'Cancer Control': 'Control'}
    )
    
    # Get available genes
    available_genes = [g for g in top_genes if g in adata_ml.var_names]
    
    if len(available_genes) < 5:
        print("  Not enough genes found in AnnData, skipping heatmap")
        return
    
    fig = plt.figure(figsize=(14, 8))
    gs = gridspec.GridSpec(1, 2, figure=fig, width_ratios=[3, 1], wspace=0.05)
    
    # Subsample for manageable heatmap
    np.random.seed(42)
    n_per_group = 200
    ctrl_idx = np.where(adata_ml.obs['condition'] == 'Control')[0]
    tb_idx = np.where(adata_ml.obs['condition'] == 'TB-affected')[0]
    
    ctrl_sample = np.random.choice(ctrl_idx, min(n_per_group, len(ctrl_idx)), replace=False)
    tb_sample = np.random.choice(tb_idx, min(n_per_group, len(tb_idx)), replace=False)
    sample_idx = np.concatenate([ctrl_sample, tb_sample])
    
    adata_sub = adata_ml[sample_idx][:, available_genes]
    
    # Get expression matrix
    if sparse.issparse(adata_sub.X):
        expr_matrix = adata_sub.X.toarray()
    else:
        expr_matrix = np.array(adata_sub.X)
    
    # Z-score normalize per gene
    from scipy.stats import zscore
    expr_z = zscore(expr_matrix, axis=0, nan_policy='omit')
    expr_z = np.nan_to_num(expr_z, 0)
    expr_z = np.clip(expr_z, -3, 3)
    
    # Main heatmap
    ax_heat = fig.add_subplot(gs[0, 0])
    im = ax_heat.imshow(expr_z.T, aspect='auto', cmap='RdBu_r', vmin=-3, vmax=3, interpolation='nearest')
    
    ax_heat.set_yticks(range(len(available_genes)))
    ax_heat.set_yticklabels(available_genes, fontsize=7, style='italic')
    ax_heat.set_xlabel('Cells (Control | TB-affected)')
    ax_heat.set_title('Expression of top 30 SHAP signature genes', fontweight='bold', fontsize=12)
    
    # Add condition color bar on top
    n_ctrl = len(ctrl_sample)
    n_tb = len(tb_sample)
    condition_colors = [CTRL_COLOR] * n_ctrl + [TB_COLOR] * n_tb
    
    # Add separator line
    ax_heat.axvline(x=n_ctrl - 0.5, color='black', linewidth=2)
    
    # Add text labels
    ax_heat.text(n_ctrl/2, -1.5, 'Control', ha='center', fontsize=9, fontweight='bold', color=CTRL_COLOR)
    ax_heat.text(n_ctrl + n_tb/2, -1.5, 'TB-affected', ha='center', fontsize=9, fontweight='bold', color=TB_COLOR)
    
    plt.colorbar(im, ax=ax_heat, label='Z-score', shrink=0.6, pad=0.02)
    
    # Direction annotation
    ax_dir = fig.add_subplot(gs[0, 1])
    sig_sub = signature[signature['gene'].isin(available_genes)].set_index('gene').loc[available_genes]
    dir_colors = [TB_COLOR if d == 'UP' else CTRL_COLOR for d in sig_sub['direction']]
    ax_dir.barh(range(len(available_genes)), sig_sub['mean_abs_shap'], color=dir_colors, 
                edgecolor='black', linewidth=0.3)
    ax_dir.set_yticks([])
    ax_dir.set_xlabel('Mean |SHAP|')
    ax_dir.set_title('Importance', fontweight='bold')
    ax_dir.invert_yaxis()
    
    legend_elements = [Patch(facecolor=TB_COLOR, label='Up in TB'),
                       Patch(facecolor=CTRL_COLOR, label='Down in TB')]
    ax_dir.legend(handles=legend_elements, loc='lower right', frameon=True, edgecolor='gray', fontsize=7)
    
    plt.savefig(os.path.join(FIG_DIR, "figure5_expression_heatmap.png"), dpi=300, facecolor='white')
    plt.savefig(os.path.join(FIG_DIR, "figure5_expression_heatmap.pdf"), facecolor='white')
    plt.close()
    print("  Saved Figure 5")


# ================================================================
# SUPPLEMENTARY: Individual high-res copies of key plots
# ================================================================
def supplementary_individual_plots(X_train, y_train, pca_model):
    """Generate individual high-resolution versions of key panels."""
    print("\n--- Generating Supplementary Figures ---")
    
    # Supp 1: PCA 3D visualization (PC1 vs PC2 vs PC3)
    from mpl_toolkits.mplot3d import Axes3D
    
    X_pca = pca_model.transform(X_train)
    
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')
    ctrl = y_train == 0
    tb = y_train == 1
    
    # Subsample for clarity
    np.random.seed(42)
    n_sub = 3000
    ctrl_idx = np.random.choice(np.where(ctrl)[0], min(n_sub, ctrl.sum()), replace=False)
    tb_idx = np.random.choice(np.where(tb)[0], min(n_sub, tb.sum()), replace=False)
    
    ax.scatter(X_pca[tb_idx, 0], X_pca[tb_idx, 1], X_pca[tb_idx, 2],
               c=TB_COLOR, alpha=0.15, s=3, label='TB-affected')
    ax.scatter(X_pca[ctrl_idx, 0], X_pca[ctrl_idx, 1], X_pca[ctrl_idx, 2],
               c=CTRL_COLOR, alpha=0.5, s=6, label='Control')
    ax.set_xlabel(f'PC1 ({pca_model.explained_variance_ratio_[0]*100:.1f}%)')
    ax.set_ylabel(f'PC2 ({pca_model.explained_variance_ratio_[1]*100:.1f}%)')
    ax.set_zlabel(f'PC3 ({pca_model.explained_variance_ratio_[2]*100:.1f}%)')
    ax.set_title('3D PCA of ML features', fontweight='bold')
    ax.legend(markerscale=3)
    
    plt.savefig(os.path.join(FIG_DIR, "supp_pca_3d.png"), dpi=300, facecolor='white')
    plt.close()
    
    # Supp 2: Correlation matrix of top 15 features
    signature = pd.read_csv(os.path.join(RESULTS_DIR, "disease_signature.csv"))
    top15 = signature.head(15)['gene'].values
    gene_names_df = pd.read_csv(os.path.join(DATA_DIR, "ml_ready", "gene_names.csv"))
    gene_list = gene_names_df.iloc[:, 0].values if gene_names_df.shape[1] == 1 else gene_names_df.iloc[:, 1].values
    
    gene_indices = []
    gene_labels = []
    for g in top15:
        if g in gene_list:
            gene_indices.append(np.where(gene_list == g)[0][0])
            gene_labels.append(g)
    
    if len(gene_indices) > 5:
        X_sub = X_train[:, gene_indices]
        corr = np.corrcoef(X_sub.T)
        
        fig, ax = plt.subplots(figsize=(8, 7))
        mask = np.triu(np.ones_like(corr, dtype=bool), k=1)
        sns.heatmap(corr, mask=mask, xticklabels=gene_labels, yticklabels=gene_labels,
                    cmap='RdBu_r', center=0, vmin=-1, vmax=1, annot=True, fmt='.2f',
                    annot_kws={'size': 7}, square=True, linewidths=0.5, ax=ax)
        ax.set_title('Pairwise correlation of top SHAP features', fontweight='bold')
        plt.xticks(rotation=45, ha='right', fontsize=8, style='italic')
        plt.yticks(fontsize=8, style='italic')
        
        plt.savefig(os.path.join(FIG_DIR, "supp_correlation_matrix.png"), dpi=300, facecolor='white')
        plt.close()
    
    print("  Saved supplementary figures")


# ================================================================
# MAIN
# ================================================================
def main():
    print("=" * 60)
    print("GENERATING MANUSCRIPT FIGURES")
    print("=" * 60)
    
    # Load all data
    adata, X_train, y_train, X_test, y_test, train_cells, test_cells, gene_names = load_data()
    
    # Generate all figures
    figure1_data_overview(adata, y_train, y_test)
    pca_model = figure2_dimensionality_reduction(adata, X_train, y_train, X_test, y_test)
    figure3_model_performance(X_train, y_train, X_test, y_test)
    figure4_shap_signature(gene_names)
    figure5_expression_heatmap(adata, gene_names)
    figure6_drug_repurposing()
    supplementary_individual_plots(X_train, y_train, pca_model)
    
    print("\n" + "=" * 60)
    print(f"ALL FIGURES SAVED TO: {FIG_DIR}")
    print("=" * 60)
    
    # List all generated files
    for f in sorted(os.listdir(FIG_DIR)):
        fpath = os.path.join(FIG_DIR, f)
        size_kb = os.path.getsize(fpath) / 1024
        print(f"  {f:45s} ({size_kb:.0f} KB)")


if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    main()
