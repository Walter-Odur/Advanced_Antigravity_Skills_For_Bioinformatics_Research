"""
Step 01: Load Raw Data & Quality Control
=========================================
Loads the Seq-Well single-cell expression matrix and metadata from the
Broad Institute SCP3227 dataset (TB lung granulomas, Sub-Saharan African cohort).

Performs standard scRNA-seq QC:
  - Filters low-quality cells/genes
  - Removes high-mitochondrial-content cells
  - Normalizes, log-transforms, and selects highly variable genes
  - Saves processed AnnData object for downstream ML

Usage:
    python 01_load_and_qc.py
"""

import scanpy as sc
import pandas as pd
import numpy as np
import os
import warnings

warnings.filterwarnings("ignore", category=FutureWarning)

# ---------- Configuration ----------
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=120, facecolor="white")
FIGDIR = "results/figures/qc"
os.makedirs(FIGDIR, exist_ok=True)
sc.settings.figdir = FIGDIR

# QC thresholds
MIN_GENES_PER_CELL = 200
MIN_CELLS_PER_GENE = 3
MAX_PCT_MITO = 20  # TB lung tissue -- higher MT threshold than PBMCs

DATA_DIR = "data"
RAW_COUNTS = os.path.join(DATA_DIR, "final_sw_dataset_raw_counts.csv")
METADATA = os.path.join(DATA_DIR, "TB_single_cell_meta_UMAP.csv")
OUTPUT_H5AD = os.path.join(DATA_DIR, "processed_tb_lung.h5ad")


def load_data():
    """Load expression matrix and metadata, handling SCP format quirks."""
    print("=" * 60)
    print("STEP 1: Loading raw data")
    print("=" * 60)

    # --- Load metadata ---
    print("\n[1/4] Loading metadata...")
    meta = pd.read_csv(METADATA, index_col=0)

    # SCP metadata has a TYPE row as the first data row (describes column types)
    # e.g., 'numeric', 'group' -- must be dropped
    if meta.index[0] == "TYPE":
        print("  -> Detected SCP 'TYPE' header row -- removing it.")
        meta = meta.iloc[1:]

    # Convert numeric columns from string to proper dtypes
    numeric_cols = ["nCount_RNA", "nFeature_RNA", "percent_mito", "percent_ribo",
                    "UMAP_1", "UMAP_2"]
    for col in numeric_cols:
        if col in meta.columns:
            meta[col] = pd.to_numeric(meta[col], errors="coerce")

    print(f"  -> Metadata: {meta.shape[0]} cells x {meta.shape[1]} annotations")

    # --- Load expression matrix ---
    print("\n[2/4] Loading raw count matrix (this may take a few minutes)...")
    counts = pd.read_csv(RAW_COUNTS, index_col=0)
    print(f"  -> Raw matrix shape: {counts.shape}")

    # Determine orientation: genes should be features (columns), cells should be rows
    # Check which axis overlaps with metadata cell barcodes
    print("\n[3/4] Orienting matrix (cells x genes)...")
    row_overlap = len(counts.index.intersection(meta.index))
    col_overlap = len(counts.columns.intersection(meta.index))
    print(f"  -> Row overlap with metadata: {row_overlap}")
    print(f"  -> Column overlap with metadata: {col_overlap}")

    if col_overlap > row_overlap:
        # Cell barcodes are in columns -- genes are rows -- need to transpose
        print("  -> Transposing: genes (rows) --> features (columns)")
        counts = counts.T
    elif row_overlap == 0 and col_overlap == 0:
        raise ValueError("No overlap between count matrix and metadata cell IDs!")
    print(f"  -> Oriented matrix: {counts.shape[0]} cells x {counts.shape[1]} genes")

    # --- Align metadata with count matrix ---
    print("\n[4/4] Aligning metadata with expression matrix...")
    common_cells = counts.index.intersection(meta.index)
    print(f"  -> {len(common_cells)} cells in common "
          f"(counts: {counts.shape[0]}, meta: {meta.shape[0]})")

    counts = counts.loc[common_cells]
    meta = meta.loc[common_cells]

    # --- Create AnnData ---
    adata = sc.AnnData(X=counts.values.astype(np.float32),
                       obs=meta,
                       var=pd.DataFrame(index=counts.columns))
    adata.obs_names = pd.Index(common_cells)
    adata.var_names = pd.Index(counts.columns)

    print(f"\n[OK] AnnData created: {adata.n_obs} cells x {adata.n_vars} genes")
    print(f"  Disease_Status distribution:\n{adata.obs['Disease_Status'].value_counts().to_string()}")
    return adata


def perform_qc(adata):
    """Run standard single-cell QC: filter, normalize, select HVGs."""
    print("\n" + "=" * 60)
    print("STEP 2: Quality Control & Preprocessing")
    print("=" * 60)

    n_before = adata.n_obs
    print(f"\n  Starting cells: {n_before}")

    # --- Cell & gene filtering ---
    sc.pp.filter_cells(adata, min_genes=MIN_GENES_PER_CELL)
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS_PER_GENE)
    print(f"  After min_genes={MIN_GENES_PER_CELL} & min_cells={MIN_CELLS_PER_GENE}: "
          f"{adata.n_obs} cells, {adata.n_vars} genes")

    # --- Mitochondrial QC ---
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"],
                               percent_top=None, log1p=False, inplace=True)

    # Save QC violin plots
    try:
        sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
                     jitter=0.4, multi_panel=True, show=False,
                     save="_qc_metrics.png")
    except Exception:
        print("  (Skipped QC violin plot -- no display available)")

    adata = adata[adata.obs.pct_counts_mt < MAX_PCT_MITO, :].copy()
    print(f"  After MT% < {MAX_PCT_MITO}%: {adata.n_obs} cells")
    print(f"  -> Removed {n_before - adata.n_obs} low-quality cells "
          f"({(n_before - adata.n_obs) / n_before * 100:.1f}%)")

    # --- Normalization ---
    print("\n  Normalizing (target_sum=10,000) and log-transforming...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # --- Highly Variable Genes ---
    print("  Identifying highly variable genes...")
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    n_hvg = adata.var.highly_variable.sum()
    print(f"  -> {n_hvg} highly variable genes identified")

    # Store raw (all genes, normalized) for later use in SHAP/DE
    adata.raw = adata

    # Subset to HVGs for downstream ML
    adata = adata[:, adata.var.highly_variable].copy()
    print(f"\n[OK] Final processed dataset: {adata.n_obs} cells x {adata.n_vars} HVGs")
    return adata


if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    os.makedirs("results/figures/qc", exist_ok=True)

    adata = load_data()
    adata_qc = perform_qc(adata)

    print(f"\nSaving processed AnnData to '{OUTPUT_H5AD}'...")
    adata_qc.write(OUTPUT_H5AD)
    print("[OK] Done! Ready for ML pipeline (Step 02).")
