"""
Step 05: SHAP Analysis -- Disease Signature Extraction
=======================================================
Uses SHAP (SHapley Additive exPlanations) to interpret the best-performing
ML model from Step 03 and extract the TB disease signature genes.

Per ml-engineer skill:
  - TreeExplainer provides exact SHAP values for tree-based models
  - Feature importance via mean |SHAP value| is more reliable than
    permutation importance or Gini importance

Outputs:
  - SHAP summary (beeswarm) plot
  - SHAP bar (feature importance) plot
  - SHAP dependence plots for top genes
  - UP/DOWN gene lists for CMap drug repurposing
  - Disease signature CSV

Usage:
    python 05_shap_analysis.py
"""

import numpy as np
import pandas as pd
import os
import json
import joblib
import warnings

import shap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# ---------- Configuration ----------
DATA_DIR = os.path.join("data", "ml_ready")
MODELS_DIR = "models"
RESULTS_DIR = "results"
FIGURES_DIR = os.path.join(RESULTS_DIR, "figures", "shap")

TOP_N_GENES = 100  # Top N genes for disease signature
TOP_N_DEPENDENCE = 5  # Top N genes for dependence plots
MAX_DISPLAY = 30  # Max genes in summary plots

# For large datasets, subsample for SHAP computation
SHAP_SAMPLE_SIZE = 2000  # Use a subsample for faster SHAP if dataset is large


def load_best_model_and_data():
    """Load the best model and test data."""
    print("=" * 60)
    print("Loading Best Model & Data")
    print("=" * 60)

    # Read which model was best
    with open(os.path.join(MODELS_DIR, "best_model.txt"), "r") as f:
        best_name = f.read().strip()
    print(f"  Best model: {best_name}")

    # Load model
    model = joblib.load(os.path.join(MODELS_DIR, f"{best_name}.joblib"))
    print(f"  Model loaded: {type(model).__name__}")

    # Load data
    train = np.load(os.path.join(DATA_DIR, "train_data.npz"))
    test = np.load(os.path.join(DATA_DIR, "test_data.npz"))
    gene_names = pd.read_csv(os.path.join(DATA_DIR, "gene_names.csv"))["gene"].tolist()

    X_train, y_train = train["X"], train["y"]
    X_test, y_test = test["X"], test["y"]

    print(f"  Train: {X_train.shape}, Test: {X_test.shape}")
    print(f"  Features: {len(gene_names)} genes")

    return model, best_name, X_train, X_test, y_train, y_test, gene_names


def compute_shap_values(model, best_name, X_train, X_test, gene_names):
    """Compute SHAP values using TreeExplainer with interventional mode.
    
    Uses interventional SHAP (Pearl's do-calculus) to handle correlated gene
    expression features correctly (Lundberg et al., Nat Mach Intell 2020;
    Janzing et al., AISTATS 2020).
    """
    print("\n" + "=" * 60)
    print("Computing SHAP Values (Interventional Mode)")
    print("=" * 60)

    # Subsample for faster SHAP computation if needed
    if X_test.shape[0] > SHAP_SAMPLE_SIZE:
        rng = np.random.RandomState(42)
        idx = rng.choice(X_test.shape[0], SHAP_SAMPLE_SIZE, replace=False)
        X_explain = X_test[idx]
        print(f"  Subsampled {SHAP_SAMPLE_SIZE} cells from test set for SHAP")
    else:
        X_explain = X_test
        print(f"  Using full test set ({X_test.shape[0]} cells) for SHAP")

    # Background dataset for interventional SHAP (500 samples — more is diminishing returns)
    BACKGROUND_SIZE = 500
    if X_train.shape[0] > BACKGROUND_SIZE:
        background = X_train[:BACKGROUND_SIZE]
    else:
        background = X_train
    print(f"  Background dataset: {background.shape[0]} cells (for interventional SHAP)")

    # For Stacking classifier, use the best base model instead
    if best_name == "StackingEnsemble":
        print("  -> StackingEnsemble detected. Using best base estimator for TreeExplainer...")
        # Try to find which base model has highest individual AUC
        # Fall back to first tree-based estimator
        base_model = None
        for name, est in model.named_estimators_.items():
            if hasattr(est, "feature_importances_"):
                base_model = est
                print(f"    Using base estimator: {name} ({type(est).__name__})")
                break

        if base_model is None:
            # Fall back to loading the best non-stacking model
            print("    Falling back to best non-stacking model...")
            test_results = pd.read_csv(os.path.join(RESULTS_DIR, "test_results.csv"),
                                       index_col=0)
            non_stacking = test_results.drop("StackingEnsemble", errors="ignore")
            fallback_name = non_stacking["roc_auc"].idxmax()
            base_model = joblib.load(os.path.join(MODELS_DIR, f"{fallback_name}.joblib"))
            print(f"    Using {fallback_name}")

        explainer = shap.TreeExplainer(
            base_model,
            data=background,
            feature_perturbation="interventional",
            model_output="probability",
        )
    else:
        explainer = shap.TreeExplainer(
            model,
            data=background,
            feature_perturbation="interventional",
            model_output="probability",
        )

    print("  Computing SHAP values (this may take a few minutes)...")
    shap_values = explainer.shap_values(X_explain)

    # For binary classification, shap_values may be a list [class0, class1]
    # We want the values for the positive class (TB-affected = class 1)
    if isinstance(shap_values, list):
        shap_values = shap_values[1]  # Class 1 = TB-affected
        print("  -> Extracted SHAP values for class 1 (TB-affected)")
    elif shap_values.ndim == 3:
        shap_values = shap_values[:, :, 1]
        print("  -> Extracted SHAP values for class 1 (TB-affected)")

    print(f"  SHAP values shape: {shap_values.shape}")

    # Cache SHAP values for reuse
    np.save(os.path.join(RESULTS_DIR, "shap_values_disease.npy"), shap_values)
    print(f"  -> Cached SHAP values to {RESULTS_DIR}/shap_values_disease.npy")

    return shap_values, X_explain


def extract_disease_signature(shap_values, gene_names):
    """Extract the TB disease signature from SHAP values."""
    print("\n" + "=" * 60)
    print("Extracting Disease Signature")
    print("=" * 60)

    # Mean absolute SHAP value per gene = global feature importance
    mean_abs_shap = np.abs(shap_values).mean(axis=0)

    # Mean SHAP value (signed) = direction of effect
    mean_shap = shap_values.mean(axis=0)

    # Create signature DataFrame
    signature = pd.DataFrame({
        "gene": gene_names,
        "mean_abs_shap": mean_abs_shap,
        "mean_shap": mean_shap,
        "direction": np.where(mean_shap > 0, "UP", "DOWN"),
    })
    signature = signature.sort_values("mean_abs_shap", ascending=False)
    signature = signature.reset_index(drop=True)
    signature["rank"] = range(1, len(signature) + 1)

    # Top N disease signature genes
    top_genes = signature.head(TOP_N_GENES)

    n_up = (top_genes["direction"] == "UP").sum()
    n_down = (top_genes["direction"] == "DOWN").sum()

    print(f"  Top {TOP_N_GENES} disease signature genes:")
    print(f"    -> {n_up} UPREGULATED in TB")
    print(f"    -> {n_down} DOWNREGULATED in TB")
    print(f"\n  Top 20 genes by |SHAP|:")
    print(top_genes[["rank", "gene", "mean_abs_shap", "direction"]].head(20).to_string(index=False))

    return signature, top_genes


def create_gene_lists(top_genes):
    """Create UP and DOWN gene lists for CMap query."""
    print("\n" + "=" * 60)
    print("Creating Gene Lists for Drug Repurposing")
    print("=" * 60)

    up_genes = top_genes[top_genes["direction"] == "UP"]["gene"].tolist()
    down_genes = top_genes[top_genes["direction"] == "DOWN"]["gene"].tolist()

    print(f"  UP genes ({len(up_genes)}):   {', '.join(up_genes[:10])}...")
    print(f"  DOWN genes ({len(down_genes)}): {', '.join(down_genes[:10])}...")

    return up_genes, down_genes


def plot_shap_results(shap_values, X_explain, gene_names, signature):
    """Generate SHAP visualization plots."""
    print("\n" + "=" * 60)
    print("Generating SHAP Plots")
    print("=" * 60)

    os.makedirs(FIGURES_DIR, exist_ok=True)

    feature_names = np.array(gene_names)

    # --- 1. SHAP Summary (Beeswarm) Plot ---
    print("  [1/4] SHAP beeswarm plot...")
    fig, ax = plt.subplots(figsize=(12, 10))
    shap.summary_plot(shap_values, X_explain, feature_names=feature_names,
                      max_display=MAX_DISPLAY, show=False)
    plt.title("SHAP Feature Importance -- TB Disease Signature\n"
              "(Positive = drives TB prediction, Negative = drives Control prediction)",
              fontsize=13, fontweight="bold", pad=15)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "shap_beeswarm.png"), dpi=200, bbox_inches="tight")
    plt.close("all")

    # --- 2. SHAP Bar Plot ---
    print("  [2/4] SHAP bar plot...")
    fig, ax = plt.subplots(figsize=(10, 8))
    shap.summary_plot(shap_values, X_explain, feature_names=feature_names,
                      plot_type="bar", max_display=MAX_DISPLAY, show=False)
    plt.title("Mean |SHAP Value| -- Gene Importance for TB Classification",
              fontsize=13, fontweight="bold", pad=15)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "shap_bar.png"), dpi=200, bbox_inches="tight")
    plt.close("all")

    # --- 3. SHAP Dependence Plots for Top Genes ---
    print(f"  [3/4] SHAP dependence plots (top {TOP_N_DEPENDENCE})...")
    top_gene_indices = signature.head(TOP_N_DEPENDENCE).index.tolist()

    fig, axes = plt.subplots(1, TOP_N_DEPENDENCE, figsize=(5 * TOP_N_DEPENDENCE, 5))
    if TOP_N_DEPENDENCE == 1:
        axes = [axes]

    for i, gene_idx in enumerate(top_gene_indices):
        gene_name = signature.loc[gene_idx, "gene"]
        col_idx = gene_names.index(gene_name)
        ax = axes[i]
        shap.dependence_plot(col_idx, shap_values, X_explain,
                             feature_names=feature_names,
                             ax=ax, show=False)
        ax.set_title(gene_name, fontsize=12, fontweight="bold")

    fig.suptitle("SHAP Dependence -- Top TB Signature Genes", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "shap_dependence.png"), dpi=200, bbox_inches="tight")
    plt.close("all")

    # --- 4. Custom Direction Plot for Top Genes ---
    print("  [4/4] Signature direction plot...")
    top_for_plot = signature.head(30).copy()
    top_for_plot = top_for_plot.sort_values("mean_shap")

    fig, ax = plt.subplots(figsize=(10, 10))
    colors = ["#2196F3" if d == "UP" else "#F44336" for d in top_for_plot["direction"]]
    ax.barh(range(len(top_for_plot)), top_for_plot["mean_shap"], color=colors,
            edgecolor="white", linewidth=0.5)
    ax.set_yticks(range(len(top_for_plot)))
    ax.set_yticklabels(top_for_plot["gene"], fontsize=9)
    ax.set_xlabel("Mean SHAP Value", fontsize=12)
    ax.set_title("TB Disease Signature -- Top 30 Genes\n"
                 "(Blue = Upregulated in TB, Red = Downregulated in TB)",
                 fontsize=13, fontweight="bold")
    ax.axvline(x=0, color="black", linewidth=0.8)
    ax.grid(axis="x", alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "signature_direction.png"), dpi=200,
                bbox_inches="tight")
    plt.close("all")

    print(f"  -> Plots saved to '{FIGURES_DIR}/'")


def save_signature(signature, top_genes, up_genes, down_genes):
    """Save disease signature and gene lists."""
    print("\n" + "=" * 60)
    print("Saving Disease Signature")
    print("=" * 60)

    os.makedirs(RESULTS_DIR, exist_ok=True)

    # Full signature (all genes ranked)
    signature.to_csv(os.path.join(RESULTS_DIR, "full_gene_ranking.csv"), index=False)

    # Top disease signature
    top_genes.to_csv(os.path.join(RESULTS_DIR, "disease_signature.csv"), index=False)

    # UP and DOWN gene lists (one gene per line, for CMap input)
    up_path = os.path.join(RESULTS_DIR, "signature_UP_genes.txt")
    down_path = os.path.join(RESULTS_DIR, "signature_DOWN_genes.txt")

    with open(up_path, "w") as f:
        f.write("\n".join(up_genes))
    with open(down_path, "w") as f:
        f.write("\n".join(down_genes))

    # Summary
    summary = {
        "n_signature_genes": len(top_genes),
        "n_up": len(up_genes),
        "n_down": len(down_genes),
        "top_10_genes": top_genes["gene"].head(10).tolist(),
    }
    with open(os.path.join(RESULTS_DIR, "shap_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    print(f"  -> disease_signature.csv ({len(top_genes)} genes)")
    print(f"  -> signature_UP_genes.txt ({len(up_genes)} genes)")
    print(f"  -> signature_DOWN_genes.txt ({len(down_genes)} genes)")
    print(f"  -> full_gene_ranking.csv (all {len(signature)} genes)")

    print(f"\n[OK] Disease signature extracted!")


if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    # Load model and data
    model, best_name, X_train, X_test, y_train, y_test, gene_names = \
        load_best_model_and_data()

    # Compute SHAP values
    shap_values, X_explain = compute_shap_values(model, best_name,
                                                  X_train, X_test, gene_names)

    # Extract disease signature
    signature, top_genes = extract_disease_signature(shap_values, gene_names)

    # Create UP/DOWN gene lists
    up_genes, down_genes = create_gene_lists(top_genes)

    # Generate plots
    plot_shap_results(shap_values, X_explain, gene_names, signature)

    # Save everything
    save_signature(signature, top_genes, up_genes, down_genes)

    print("\n--> Next: Run 06_drug_repurposing.py to find candidate drugs.")
