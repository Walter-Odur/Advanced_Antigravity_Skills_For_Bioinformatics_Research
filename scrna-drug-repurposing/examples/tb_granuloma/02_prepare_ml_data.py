"""
Step 02: Prepare ML-Ready Data
================================
Takes the QC'd AnnData from Step 01 and creates the ML classification dataset.

Key decisions (guided by ML-engineer & scikit-learn best practices):
  - Binary label: TB-affected (TB + HIVTB) vs. Control (Cancer Control)
  - HIV-only cells excluded (confounding without TB)
  - Stratified train/test split (80/20) preserving class ratios
  - Feature matrix = HVG-filtered log-normalized counts

Usage:
    python 02_prepare_ml_data.py
"""

import scanpy as sc
import pandas as pd
import numpy as np
import os
import joblib
from scipy import sparse
from sklearn.model_selection import train_test_split

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

# ---------- Configuration ----------
DATA_DIR = "data"
INPUT_H5AD = os.path.join(DATA_DIR, "processed_tb_lung.h5ad")
OUTPUT_DIR = os.path.join(DATA_DIR, "ml_ready")

# Label mapping: TB-positive vs. Control
# HIVTB and TB --> 1 (TB-affected)
# Cancer Control --> 0 (Non-TB control)
# HIV-only --> EXCLUDED (confounds TB signal without active TB)
LABEL_MAP = {
    "TB": 1,
    "HIVTB": 1,
    "Cancer Control": 0,
}
EXCLUDED_GROUPS = ["HIV"]

TEST_SIZE = 0.20
RANDOM_STATE = 42


def prepare_labels(adata):
    """Create binary TB classification labels and filter excluded groups."""
    print("=" * 60)
    print("STEP 1: Preparing Classification Labels")
    print("=" * 60)

    disease_col = "Disease_Status"
    print(f"\n  Original Disease_Status distribution:")
    print(f"  {adata.obs[disease_col].value_counts().to_string()}")

    # Exclude confounding groups
    mask_exclude = adata.obs[disease_col].isin(EXCLUDED_GROUPS)
    n_excluded = mask_exclude.sum()
    print(f"\n  Excluding {n_excluded} cells from groups: {EXCLUDED_GROUPS}")
    adata = adata[~mask_exclude].copy()

    # Map to binary labels
    adata.obs["label"] = adata.obs[disease_col].map(LABEL_MAP)

    # Sanity check: no NaN labels
    n_unmapped = adata.obs["label"].isna().sum()
    if n_unmapped > 0:
        unmapped = adata.obs.loc[adata.obs["label"].isna(), disease_col].unique()
        raise ValueError(f"Unmapped Disease_Status values: {unmapped}")

    adata.obs["label"] = adata.obs["label"].astype(int)

    print(f"\n  Binary label distribution:")
    label_counts = adata.obs["label"].value_counts()
    for label, count in label_counts.items():
        name = "TB-affected" if label == 1 else "Control"
        print(f"    {name} (label={label}): {count} cells")

    ratio = label_counts.max() / label_counts.min()
    print(f"  -> Class imbalance ratio: {ratio:.1f}:1")
    print(f"    (Will use class_weight='balanced' & stratified CV to handle this)")

    return adata


def create_feature_matrix(adata):
    """Extract the feature matrix (X) from AnnData."""
    print("\n" + "=" * 60)
    print("STEP 2: Creating Feature Matrix")
    print("=" * 60)

    # Extract the dense matrix (HVG-filtered, log-normalized)
    if sparse.issparse(adata.X):
        X = np.asarray(adata.X.todense(), dtype=np.float32)
    else:
        X = np.asarray(adata.X, dtype=np.float32)

    y = adata.obs["label"].values
    gene_names = adata.var_names.tolist()
    cell_ids = adata.obs_names.tolist()

    print(f"  Feature matrix X: {X.shape[0]} cells x {X.shape[1]} genes")
    print(f"  Label vector y:   {len(y)} labels")
    print(f"  Gene features:    {len(gene_names)} HVGs")

    return X, y, gene_names, cell_ids


def split_data(X, y, cell_ids):
    """Stratified train/test split preserving class ratios."""
    print("\n" + "=" * 60)
    print("STEP 3: Train/Test Split")
    print("=" * 60)

    X_train, X_test, y_train, y_test, idx_train, idx_test = train_test_split(
        X, y, np.arange(len(y)),
        test_size=TEST_SIZE,
        stratify=y,
        random_state=RANDOM_STATE,
    )

    train_cells = [cell_ids[i] for i in idx_train]
    test_cells = [cell_ids[i] for i in idx_test]

    print(f"  Train set: {X_train.shape[0]} cells ({X_train.shape[0]/len(y)*100:.1f}%)")
    print(f"    -> label 0: {(y_train == 0).sum()}, label 1: {(y_train == 1).sum()}")
    print(f"  Test set:  {X_test.shape[0]} cells ({X_test.shape[0]/len(y)*100:.1f}%)")
    print(f"    -> label 0: {(y_test == 0).sum()}, label 1: {(y_test == 1).sum()}")

    return X_train, X_test, y_train, y_test, train_cells, test_cells


def save_ml_data(X_train, X_test, y_train, y_test, gene_names,
                 train_cells, test_cells):
    """Save ML-ready arrays and metadata."""
    print("\n" + "=" * 60)
    print("STEP 4: Saving ML-Ready Data")
    print("=" * 60)

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Save as compressed numpy arrays for fast loading
    np.savez_compressed(os.path.join(OUTPUT_DIR, "train_data.npz"),
                        X=X_train, y=y_train)
    np.savez_compressed(os.path.join(OUTPUT_DIR, "test_data.npz"),
                        X=X_test, y=y_test)

    # Save gene names (feature names for SHAP)
    pd.Series(gene_names, name="gene").to_csv(
        os.path.join(OUTPUT_DIR, "gene_names.csv"), index=False)

    # Save cell IDs for traceability
    pd.Series(train_cells, name="cell_id").to_csv(
        os.path.join(OUTPUT_DIR, "train_cells.csv"), index=False)
    pd.Series(test_cells, name="cell_id").to_csv(
        os.path.join(OUTPUT_DIR, "test_cells.csv"), index=False)

    # Save a summary
    summary = {
        "n_train": len(y_train),
        "n_test": len(y_test),
        "n_features": len(gene_names),
        "train_pos": int((y_train == 1).sum()),
        "train_neg": int((y_train == 0).sum()),
        "test_pos": int((y_test == 1).sum()),
        "test_neg": int((y_test == 0).sum()),
    }
    pd.Series(summary).to_json(os.path.join(OUTPUT_DIR, "data_summary.json"))

    sizes = {}
    for f in os.listdir(OUTPUT_DIR):
        fpath = os.path.join(OUTPUT_DIR, f)
        sizes[f] = f"{os.path.getsize(fpath) / 1e6:.1f} MB"

    print(f"  Saved to '{OUTPUT_DIR}/':")
    for fname, size in sizes.items():
        print(f"    {fname}: {size}")

    print(f"\n[OK] ML data ready! ({len(y_train) + len(y_test)} total cells, "
          f"{len(gene_names)} features)")


if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    # Load processed AnnData
    print(f"Loading processed AnnData from '{INPUT_H5AD}'...\n")
    adata = sc.read_h5ad(INPUT_H5AD)
    print(f"  Loaded: {adata.n_obs} cells x {adata.n_vars} genes\n")

    # Pipeline
    adata = prepare_labels(adata)
    X, y, gene_names, cell_ids = create_feature_matrix(adata)
    X_train, X_test, y_train, y_test, train_cells, test_cells = split_data(X, y, cell_ids)
    save_ml_data(X_train, X_test, y_train, y_test, gene_names, train_cells, test_cells)

    print("\n--> Next: Run 03_ensemble_ml.py to train the ML ensemble.")
