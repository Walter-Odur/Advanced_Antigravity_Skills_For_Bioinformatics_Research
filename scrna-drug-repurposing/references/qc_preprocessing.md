# scRNA-seq QC, Preprocessing & Scanpy — Comprehensive Guide

> [!NOTE]
> This guide covers the complete Scanpy workflow for scRNA-seq analysis from raw data to ML-ready features. All methodological choices grounded in high-impact peer-reviewed literature from *Nature Methods*, *Nature Reviews Genetics*, *Genome Biology*, and *Molecular Systems Biology*.

---

## 1. Key References

1. **Luecken & Theis (2019)** — "Current best practices in scRNA-seq analysis: a tutorial." *Molecular Systems Biology*, 15(6):e8746
2. **Heumos et al. (2023)** — "Best practices for single-cell analysis across modalities." *Nature Reviews Genetics*, 24:550–572
3. **Fleming et al. (2023)** — "Unsupervised removal of systematic background noise using CellBender." *Nature Methods*, 20:1323–1335
4. **Young & Behjati (2020)** — "SoupX removes ambient RNA contamination." *GigaScience*, 9(12):giaa151
5. **Germain et al. (2021)** — "Doublet identification using scDblFinder." *F1000Research*, 10:979
6. **Wolock et al. (2019)** — "Scrublet: computational identification of cell doublets." *Cell Systems*, 8(4):281–291
7. **Lun, Bach & Marioni (2016)** — "Pooling across cells to normalize scRNA-seq data." *Genome Biology*, 17:75
8. **Hafemeister & Satija (2019)** — "Normalization and variance stabilization using regularized NB regression." *Genome Biology*, 20:296
9. **Choudhary & Satija (2022)** — "Comparison of statistical error models for scRNA-seq." *Genome Biology*, 23:27
10. **Korsunsky et al. (2019)** — "Fast, sensitive and accurate integration with Harmony." *Nature Methods*, 16:1289–1296
11. **Lopez et al. (2018)** — "Deep generative modeling for single-cell transcriptomics." *Nature Methods*, 15:1053–1058
12. **Townes et al. (2019)** — "Feature selection for scRNA-seq based on a multinomial model." *Genome Biology*, 20:295
13. **Wolf et al. (2018)** — "SCANPY: large-scale single-cell gene expression data analysis." *Genome Biology*, 19:15
14. **Traag et al. (2019)** — "From Louvain to Leiden: guaranteeing well-connected communities." *Scientific Reports*, 9:5233
15. **McInnes et al. (2018)** — "UMAP: Uniform Manifold Approximation and Projection." *arXiv:1802.03426*

---

## 2. AnnData Structure — The Core Data Object

**Wolf et al., Genome Biology 2018**

Every Scanpy operation revolves around `AnnData`:

```python
import scanpy as sc
import pandas as pd
import numpy as np

# Core attributes
adata.X           # Expression matrix (cells × genes) — sparse or dense
adata.obs         # Cell metadata (DataFrame, n_obs × n_fields)
adata.var         # Gene metadata (DataFrame, n_vars × n_fields)
adata.uns         # Unstructured annotations (dict) — colors, parameters
adata.obsm        # Multi-dimensional cell embeddings (PCA, UMAP, t-SNE)
adata.varm        # Multi-dimensional gene embeddings (PCA loadings)
adata.obsp        # Pairwise cell relationships (neighbor graph)
adata.layers      # Alternative expression matrices (raw counts, normalized)
adata.raw         # Frozen raw data backup (immutable after assignment)

# Access names
adata.obs_names   # Cell barcodes (Index)
adata.var_names   # Gene names (Index)
adata.shape       # (n_cells, n_genes)
```

### AnnData Best Practices

```python
# Always save raw before subsetting/normalizing
adata.raw = adata.copy()

# Use layers for multiple expression representations
adata.layers['counts'] = adata.X.copy()     # Raw counts
# ... normalize adata.X ...
adata.layers['normalized'] = adata.X.copy()  # Normalized

# Subset safely (always .copy() to avoid views)
adata_subset = adata[adata.obs['cell_type'] == 'macrophage'].copy()
```

---

## 3. Data Loading

### From 10X Genomics
```python
# CellRanger output (filtered)
adata = sc.read_10x_mtx('path/to/filtered_feature_bc_matrix/')
adata = sc.read_10x_h5('path/to/filtered_feature_bc_matrix.h5')

# Multiple samples
adatas = {}
for sample in ['sample1', 'sample2', 'sample3']:
    adatas[sample] = sc.read_10x_mtx(f'data/{sample}/filtered_feature_bc_matrix/')
    adatas[sample].obs['sample'] = sample

# Concatenate
adata = sc.concat(adatas, label='sample', join='outer')
```

### From h5ad / CSV / Other
```python
adata = sc.read_h5ad('data/processed.h5ad')
adata = sc.read_csv('data/expression_matrix.csv')

# From pandas DataFrame
adata = sc.AnnData(X=df_counts.values,
                    obs=pd.DataFrame(index=df_counts.index),
                    var=pd.DataFrame(index=df_counts.columns))
```

### From Broad Single Cell Portal (SCP)
```python
meta = pd.read_csv("metadata.csv", index_col=0)
counts = pd.read_csv("expression.csv", index_col=0)

# SCP quirk: TYPE row as first data row
if meta.index[0] == "TYPE":
    meta = meta.iloc[1:]

# Auto-detect matrix orientation
row_overlap = len(counts.index.intersection(meta.index))
col_overlap = len(counts.columns.intersection(meta.index))
if col_overlap > row_overlap:
    counts = counts.T  # genes in rows → cells in rows

# Build AnnData
common_cells = counts.index.intersection(meta.index)
adata = sc.AnnData(X=counts.loc[common_cells].values.astype(np.float32),
                    obs=meta.loc[common_cells],
                    var=pd.DataFrame(index=counts.columns))
```

---

## 4. Quality Control — Tiered Pipeline

### Tier 0: Basic Gene/Cell Filtering (Always)

```python
# Configure settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')

# Identify mitochondrial, ribosomal, hemoglobin genes
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
adata.var['hb'] = adata.var_names.str.match('^HB[^(P)]')

# Calculate QC metrics
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo', 'hb'],
                            percent_top=None, log1p=False, inplace=True)

# Visualize QC distributions
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True, show=False)
plt.savefig("figures/qc/qc_violins_pre_filter.png", dpi=200, bbox_inches='tight')
plt.close()

# Joint QC plot (genes vs counts, colored by MT%)
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts',
              color='pct_counts_mt', show=False)
plt.savefig("figures/qc/qc_scatter.png", dpi=200, bbox_inches='tight')
plt.close()

# Apply filters
print(f"Cells before filtering: {adata.n_obs}")
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs.pct_counts_mt < MAX_PCT_MITO].copy()
print(f"Cells after filtering: {adata.n_obs}")
```

### Tissue-Specific MT% Thresholds

| Tissue | MAX_MT% | Reference |
|---|---|---|
| PBMCs / Blood | 5–10% | Heumos et al., *Nat Rev Genet* 2023 |
| Lung | 15–25% | Travaglini et al., *Nature* 2020 |
| Liver | 20–30% | Aizarani et al., *Nature* 2019 |
| Brain | 5–10% | Lake et al., *Nat Biotechnol* 2018 |
| Tumor (solid) | 10–20% | Tirosh et al., *Science* 2016 |
| Intestine | 20–30% | Smillie et al., *Cell* 2019 |

### Tier 1: Ambient RNA Removal (Recommended for publication)

**Fleming et al., Nature Methods 2023; Young & Behjati, GigaScience 2020**

Ambient RNA from lysed cells contaminates droplets, causing false-positive gene detection.

```python
# === CellBender (deep generative model — most robust) ===
# Run as standalone command BEFORE loading into Scanpy:
# cellbender remove-background \
#   --input raw_feature_bc_matrix.h5 \
#   --output cellbender_output.h5 \
#   --expected-cells 5000 \
#   --total-droplets-included 25000 \
#   --fpr 0.01 \
#   --epochs 150

# Then load the decontaminated matrix:
adata = sc.read_10x_h5("cellbender_output_filtered.h5")

# === SoupX (faster, lighter — R package) ===
# library(SoupX)
# sc <- load10X("path/to/cellranger/output/")
# sc <- autoEstCont(sc)
# out <- adjustCounts(sc)
# Use rpy2 to bridge R ↔ Python if needed
```

**When to use:** Always for droplet-based protocols (10X Chromium, Drop-seq). Less critical for plate-based (Smart-seq2).

### Tier 2: Doublet Detection (Recommended for publication)

**Germain et al., F1000Research 2021; Wolock et al., Cell Systems 2019**

Expected doublet rate: ~0.8% per 1,000 cells captured (10X Genomics guideline).

```python
# === Scrublet (Python-native, fast) ===
import scrublet as scr

# Expected rate based on cell loading
n_cells = adata.n_obs
expected_rate = n_cells / 1000 * 0.008  # ~0.8% per 1000 cells

scrub = scr.Scrublet(adata.X, expected_doublet_rate=min(expected_rate, 0.1))
doublet_scores, predicted_doublets = scrub.scrub_doublets(
    min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30
)
adata.obs['doublet_score'] = doublet_scores
adata.obs['predicted_doublet'] = predicted_doublets

print(f"Detected {predicted_doublets.sum()} doublets ({predicted_doublets.mean()*100:.1f}%)")
adata = adata[~adata.obs['predicted_doublet']].copy()

# === scDblFinder (best benchmarked — R via Bioconductor) ===
# library(scDblFinder)
# sce <- scDblFinder(sce)
# doublets <- colData(sce)$scDblFinder.class == "doublet"

# === SOLO (scvi-tools, deep learning) ===
# import scvi
# scvi.model.SCVI.setup_anndata(adata)
# vae = scvi.model.SCVI(adata)
# vae.train()
# solo = scvi.external.SOLO.from_scvi_model(vae)
# solo.train()
# doublet_preds = solo.predict()
```

---

## 5. Normalization

### Standard: Library-Size + Log (Default for ML pipelines)

```python
# Save raw counts first
adata.layers['counts'] = adata.X.copy()

# Normalize to 10,000 counts per cell
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Save raw (all genes) before HVG subsetting
adata.raw = adata
```

**For ML classification pipelines:** Library-size + log is recommended because tree-based models (RF, XGBoost, LightGBM) are invariant to monotonic transformations.

### Advanced: Scran Pooling-Deconvolution

**Lun, Bach & Marioni, Genome Biology 2016**

More robust for datasets with many zeros. Pools cells to estimate size factors.

```python
# Via rpy2 (R interop)
import anndata2ri
from rpy2.robjects import r
anndata2ri.activate()

r('library(scran)')
r('library(BiocParallel)')

adata_raw = adata.copy()
r.assign('sce', anndata2ri.py2rpy(adata_raw))
r('sce <- computeSumFactors(sce, clusters=quickCluster(sce))')
r('sce <- logNormCounts(sce)')
adata_scran = anndata2ri.rpy2py(r('sce'))
```

### Advanced: SCTransform v2 (Pearson Residuals)

**Choudhary & Satija, Genome Biology 2022**

```python
# Via Scanpy experimental module (Scanpy >= 1.9)
import scanpy.experimental as sce
sce.pp.normalize_pearson_residuals(adata, n_top_genes=3000)
```

### Normalization Decision Guide

| Method | Best For | Limitations |
|---|---|---|
| Library-size + log | ML pipelines, most analyses | Assumes equal library complexity |
| Scran pooling | Many zeros, rare cell types | Requires R; computationally heavier |
| SCTransform v2 | Variance stabilization, integration | Pearson residuals may not suit all ML models |

---

## 6. Feature Selection

### Standard: Highly Variable Genes (HVGs)

```python
sc.pp.highly_variable_genes(
    adata,
    min_mean=0.0125,    # Minimum mean expression
    max_mean=3,         # Maximum mean expression
    min_disp=0.5,       # Minimum dispersion
    # OR: n_top_genes=2000  # Select top N by dispersion
)

# Visualize
sc.pl.highly_variable_genes(adata, show=False)
plt.savefig("figures/qc/hvg_selection.png", dpi=200, bbox_inches='tight')
plt.close()

# Subset to HVGs (keep .raw for full gene access)
adata = adata[:, adata.var.highly_variable].copy()
print(f"Selected {adata.n_vars} HVGs")
```

### Advanced: Deviance-Based Feature Selection

**Townes et al., Genome Biology 2019**

Uses deviance from a null multinomial model. More principled than dispersion-based selection.

```python
# Available via scry (R) or manual computation
# Higher deviance = more informative gene
from scipy.stats import chi2

# Simplified deviance computation
total_counts = np.array(adata.X.sum(axis=1)).flatten()
gene_counts = np.array(adata.X.sum(axis=0)).flatten()
expected_freq = gene_counts / gene_counts.sum()

# Compute deviance per gene (multinomial null model)
# Full implementation: see Townes et al. (2019) supplementary code
```

### Batch-Aware HVG Selection

For multi-batch studies, compute HVGs per batch and take the intersection:

```python
sc.pp.highly_variable_genes(adata, n_top_genes=3000, batch_key='batch')
# Genes marked highly_variable if HVG in multiple batches
```

---

## 7. Batch Correction (Multi-Sample Studies)

### Harmony (Fast, PCA-space)

**Korsunsky et al., Nature Methods 2019**

```python
import scanpy.external as sce

# Run PCA first
sc.tl.pca(adata, svd_solver='arpack', n_comps=50)

# Harmony integration
sce.pp.harmony_integrate(adata, key='batch', basis='X_pca',
                          adjusted_basis='X_pca_harmony')

# Use corrected embeddings for downstream
sc.pp.neighbors(adata, use_rep='X_pca_harmony')
```

### scVI (Deep Generative Model — Most Robust)

**Lopez et al., Nature Methods 2018**

```python
import scvi

scvi.model.SCVI.setup_anndata(adata, batch_key="batch")
model = scvi.model.SCVI(adata, n_latent=30, n_layers=2)
model.train(max_epochs=200)

# Latent representation (batch-corrected)
adata.obsm["X_scVI"] = model.get_latent_representation()

# Use for downstream
sc.pp.neighbors(adata, use_rep='X_scVI')
```

### ComBat (Simple Linear Correction)

```python
sc.pp.combat(adata, key='batch')
```

### Batch Correction Decision Guide

| Scenario | Recommendation |
|---|---|
| Single sample, single batch | Skip |
| Multiple samples, same protocol | Harmony (fast, effective) |
| Complex multi-site study | scVI (handles confounding) |
| **For ML classification** | Apply BEFORE feature extraction if batches confound with labels |

> [!WARNING]
> If disease/control labels are perfectly confounded with batch (e.g., all disease from batch A, all control from batch B), no batch correction can fix this. The experimental design must ensure labels are mixed across batches.

---

## 8. Dimensionality Reduction

### PCA

```python
# Regress out unwanted variation (optional — assess impact first)
# sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])

# Scale to unit variance
sc.pp.scale(adata, max_value=10)

# PCA
sc.tl.pca(adata, svd_solver='arpack', n_comps=50)

# Elbow plot to determine optimal n_pcs
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, show=False)
plt.savefig("figures/qc/pca_variance_ratio.png", dpi=200, bbox_inches='tight')
plt.close()

# Typically use 20-50 PCs (where elbow flattens)
N_PCS = 30  # Adjust based on elbow plot
```

### Neighborhood Graph

```python
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=N_PCS)
```

### UMAP

**McInnes et al., arXiv 2018**

```python
sc.tl.umap(adata, min_dist=0.3, spread=1.0)

# Visualize by condition
sc.pl.umap(adata, color=['condition', 'cell_type'], show=False)
plt.savefig("figures/qc/umap_overview.png", dpi=200, bbox_inches='tight')
plt.close()
```

### t-SNE

```python
sc.tl.tsne(adata, n_pcs=N_PCS, perplexity=30)
sc.pl.tsne(adata, color='condition', show=False)
```

---

## 9. Clustering

### Leiden (Recommended)

**Traag et al., Scientific Reports 2019**

```python
# Leiden clustering (preferred over Louvain — guarantees connected communities)
sc.tl.leiden(adata, resolution=0.5)

# Visualize
sc.pl.umap(adata, color='leiden', legend_loc='on data', show=False)
plt.savefig("figures/qc/leiden_clusters.png", dpi=200, bbox_inches='tight')
plt.close()

# Try multiple resolutions to find optimal granularity
for res in [0.3, 0.5, 0.8, 1.0, 1.5]:
    sc.tl.leiden(adata, resolution=res, key_added=f'leiden_{res}')
    n_clusters = adata.obs[f'leiden_{res}'].nunique()
    print(f"  Resolution {res}: {n_clusters} clusters")
```

---

## 10. Marker Gene Identification

```python
# Find marker genes for each cluster (Wilcoxon rank-sum — recommended)
sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon',
                         use_raw=True)

# Visualize
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, show=False)
plt.savefig("figures/qc/marker_genes.png", dpi=200, bbox_inches='tight')
plt.close()

# Dot plot of top markers
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, show=False)
plt.savefig("figures/qc/marker_dotplot.png", dpi=200, bbox_inches='tight')
plt.close()

# Get results as DataFrame
markers = sc.get.rank_genes_groups_df(adata, group='0')
markers = markers[markers['pvals_adj'] < 0.05]
print(f"Cluster 0: {len(markers)} significant markers")
```

---

## 11. Cell Type Annotation

### Manual (Marker-Based)

```python
# Define known markers
marker_genes = {
    'T cells': ['CD3D', 'CD3E', 'CD3G', 'IL7R'],
    'CD8 T cells': ['CD8A', 'CD8B', 'GZMB', 'PRF1'],
    'NK cells': ['NKG7', 'GNLY', 'KLRD1'],
    'B cells': ['MS4A1', 'CD79A', 'CD19'],
    'Monocytes': ['CD14', 'LYZ', 'S100A12'],
    'Macrophages': ['CD68', 'CD163', 'MRC1'],
    'Dendritic cells': ['CLEC9A', 'CD1C', 'FCER1A'],
    'Fibroblasts': ['COL1A1', 'COL3A1', 'DCN'],
    'Epithelial': ['EPCAM', 'KRT18', 'KRT19'],
}

# Score each cell type
for ct, genes in marker_genes.items():
    genes_present = [g for g in genes if g in adata.raw.var_names]
    if genes_present:
        sc.tl.score_genes(adata, genes_present, score_name=f'{ct}_score', use_raw=True)

# Visualize marker expression on UMAP
sc.pl.umap(adata, color=list(marker_genes.keys())[:4], use_raw=True, show=False)
plt.savefig("figures/qc/marker_expression.png", dpi=200, bbox_inches='tight')
plt.close()

# Manual annotation mapping
cluster_to_celltype = {
    '0': 'CD14+ Monocytes',
    '1': 'CD4 T cells',
    '2': 'NK cells',
    # ... map all clusters
}
adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_to_celltype)
```

### Automated (Reference-Based)

```python
# Using CellTypist (automated annotation)
# import celltypist
# model = celltypist.models.download_models(model='Immune_All_Low.pkl')
# predictions = celltypist.annotate(adata, model='Immune_All_Low.pkl')
# adata = predictions.to_adata()
```

---

## 12. Differential Expression Between Conditions

```python
# Compare disease vs. control within a cell type
adata_mono = adata[adata.obs['cell_type'] == 'CD14+ Monocytes'].copy()

sc.tl.rank_genes_groups(adata_mono, groupby='condition',
                         groups=['disease'], reference='control',
                         method='wilcoxon', use_raw=True)

# Results
de_results = sc.get.rank_genes_groups_df(adata_mono, group='disease')
de_results = de_results[de_results['pvals_adj'] < 0.05]
print(f"Significant DE genes: {len(de_results)}")
```

---

## 13. Trajectory Inference

### PAGA (Partition-Based Graph Abstraction)

```python
sc.tl.paga(adata, groups='leiden')
sc.pl.paga(adata, color='leiden', show=False)
plt.savefig("figures/qc/paga.png", dpi=200, bbox_inches='tight')
plt.close()

# Initialize UMAP with PAGA positions (better layout)
sc.tl.umap(adata, init_pos='paga')
```

### Diffusion Pseudotime

```python
# Set root cell (e.g., stem/progenitor cluster)
adata.uns['iroot'] = np.flatnonzero(adata.obs['leiden'] == '0')[0]
sc.tl.diffmap(adata)
sc.tl.dpt(adata)

sc.pl.umap(adata, color='dpt_pseudotime', show=False)
plt.savefig("figures/qc/pseudotime.png", dpi=200, bbox_inches='tight')
plt.close()
```

---

## 14. Gene Set Scoring

```python
# Score cells for pathway activity
inflammatory_genes = ['IL1B', 'IL6', 'TNF', 'CXCL8', 'S100A12']
sc.tl.score_genes(adata, inflammatory_genes, score_name='inflammatory_score',
                   use_raw=True)

sc.pl.umap(adata, color='inflammatory_score', show=False)
plt.savefig("figures/qc/inflammatory_score.png", dpi=200, bbox_inches='tight')
plt.close()
```

---

## 15. Label Mapping for Binary Classification

```python
# Map disease groups to binary labels
LABEL_MAP = {
    "TB": 1, "HIVTB": 1,         # Disease positive
    "Cancer Control": 0,           # Control
}
EXCLUDED_GROUPS = ["HIV"]          # Confounders — exclude entirely

adata.obs["label"] = adata.obs["Disease_Status"].map(LABEL_MAP)

# Exclude confounding groups
adata = adata[~adata.obs["Disease_Status"].isin(EXCLUDED_GROUPS)].copy()

# Drop unmapped cells
adata = adata[adata.obs["label"].notna()].copy()
adata.obs["label"] = adata.obs["label"].astype(int)

# Sanity check
print(f"Label distribution:")
print(adata.obs["label"].value_counts())
assert adata.obs["label"].isna().sum() == 0, "Unmapped labels found!"
```

---

## 16. ML Feature Extraction

```python
import numpy as np
from sklearn.model_selection import train_test_split

# Extract dense feature matrix from HVG-filtered AnnData
if hasattr(adata.X, 'toarray'):
    X = adata.X.toarray().astype(np.float32)
else:
    X = np.array(adata.X, dtype=np.float32)

y = adata.obs["label"].values.astype(int)
gene_names = adata.var_names.tolist()

# Stratified train/test split
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.20, random_state=42, stratify=y
)

# Save
np.savez_compressed("data/ml_ready/train_data.npz", X=X_train, y=y_train)
np.savez_compressed("data/ml_ready/test_data.npz", X=X_test, y=y_test)
pd.DataFrame({"gene": gene_names}).to_csv("data/ml_ready/gene_names.csv", index=False)
```

---

## 17. Publication-Quality Scanpy Plots

```python
# High-quality defaults
sc.settings.set_figure_params(dpi=300, frameon=False, figsize=(5, 5))
sc.settings.file_format_figs = 'pdf'

# UMAP with custom styling
sc.pl.umap(adata, color='cell_type',
           palette='Set2',
           legend_loc='on data',
           legend_fontsize=12,
           legend_fontoutline=2,
           frameon=False,
           save='_publication.pdf')

# Heatmap of marker genes
sc.pl.heatmap(adata, var_names=marker_genes, groupby='cell_type',
              swap_axes=True, show_gene_labels=True,
              save='_markers.pdf')

# Dot plot
sc.pl.dotplot(adata, var_names=marker_genes, groupby='cell_type',
              save='_dotplot.pdf')

# Stacked violin
sc.pl.stacked_violin(adata, var_names=['CD3D', 'CD14', 'NKG7', 'MS4A1'],
                      groupby='cell_type', save='_stacked_violin.pdf')
```

---

## 18. Saving & Loading Results

```python
# Save processed data
adata.write('data/processed_data.h5ad')

# Export metadata
adata.obs.to_csv('results/cell_metadata.csv')
adata.var.to_csv('results/gene_metadata.csv')

# Load
adata = sc.read_h5ad('data/processed_data.h5ad')
```

---

## 19. Common Pitfalls & Best Practices

| Pitfall | Solution |
|---|---|
| Not saving raw counts | Always `adata.raw = adata` before filtering genes |
| Overly strict MT% | Use tissue-specific thresholds (see table) |
| Using Louvain over Leiden | Leiden guarantees connected communities (Traag et al. 2019) |
| Single clustering resolution | Try 0.3–1.5 range; validate with marker genes |
| Ignoring batch effects | Apply Harmony/scVI if multiple batches present |
| Not checking QC plots | Always visualize distributions before filtering |
| Wrong matrix orientation | Auto-detect by checking barcode overlap with metadata |
| Memory issues with large data | Use backed mode: `sc.read_h5ad(f, backed='r')` |
| Not saving intermediate results | Write `.h5ad` at each checkpoint |

---

## 20. Key References

1. Wolf F, et al. (2018). SCANPY. *Genome Biol*, 19:15
2. Luecken MD, Theis FJ. (2019). Best practices in scRNA-seq analysis. *Mol Syst Biol*, 15(6):e8746
3. Heumos L, et al. (2023). Best practices for single-cell analysis. *Nat Rev Genet*, 24:550–572
4. Fleming SJ, et al. (2023). CellBender. *Nat Methods*, 20:1323–1335
5. Young MD, Behjati S. (2020). SoupX. *GigaScience*, 9(12):giaa151
6. Germain P-L, et al. (2021). scDblFinder. *F1000Research*, 10:979
7. Wolock SL, et al. (2019). Scrublet. *Cell Systems*, 8(4):281–291
8. Lun ATL, et al. (2016). Scran normalization. *Genome Biol*, 17:75
9. Hafemeister C, Satija R. (2019). SCTransform. *Genome Biol*, 20:296
10. Choudhary S, Satija R. (2022). SCTransform v2. *Genome Biol*, 23:27
11. Korsunsky I, et al. (2019). Harmony. *Nat Methods*, 16:1289–1296
12. Lopez R, et al. (2018). scVI. *Nat Methods*, 15:1053–1058
13. Townes FW, et al. (2019). Feature selection for scRNA-seq. *Genome Biol*, 20:295
14. Traag VA, et al. (2019). Leiden algorithm. *Sci Rep*, 9:5233
15. McInnes L, et al. (2018). UMAP. *arXiv:1802.03426*
