"""
Step 03: Ensemble Machine Learning Training & Evaluation
==========================================================
Trains and cross-validates an ensemble of 4 ML classifiers on the
TB vs. Control single-cell classification task:

  1. Random Forest        (Bagging)
  2. XGBoost              (Gradient Boosting)
  3. LightGBM             (Gradient Boosting - alternative)
  4. Stacking Classifier  (Meta-learner combining 1-3)

Best practices applied (from ml-engineer & scikit-learn skills):
  - Stratified K-Fold cross-validation
  - class_weight='balanced' / scale_pos_weight for imbalanced classes
  - Comprehensive metrics: Accuracy, F1, Precision, Recall, AUC-ROC
  - Confusion matrices & ROC curves saved to results/
  - Models serialized via joblib

Usage:
    python 03_ensemble_ml.py
"""

import numpy as np
import pandas as pd
import os
import json
import time
import joblib
import warnings

from sklearn.ensemble import (
    RandomForestClassifier,
    StackingClassifier,
)
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_validate
from sklearn.metrics import (
    accuracy_score,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
    classification_report,
    confusion_matrix,
    RocCurveDisplay,
)
from xgboost import XGBClassifier
from lightgbm import LGBMClassifier

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# ---------- Configuration ----------
DATA_DIR = os.path.join("data", "ml_ready")
RESULTS_DIR = "results"
FIGURES_DIR = os.path.join(RESULTS_DIR, "figures", "ml")
MODELS_DIR = "models"

N_CV_FOLDS = 5
RANDOM_STATE = 42
N_JOBS = -1  # Use all CPU cores


def load_data():
    """Load the ML-ready train/test data."""
    print("=" * 60)
    print("Loading ML-Ready Data")
    print("=" * 60)

    train = np.load(os.path.join(DATA_DIR, "train_data.npz"))
    test = np.load(os.path.join(DATA_DIR, "test_data.npz"))
    gene_names = pd.read_csv(os.path.join(DATA_DIR, "gene_names.csv"))["gene"].tolist()

    X_train, y_train = train["X"], train["y"]
    X_test, y_test = test["X"], test["y"]

    print(f"  Train: {X_train.shape[0]} cells x {X_train.shape[1]} features")
    print(f"  Test:  {X_test.shape[0]} cells x {X_test.shape[1]} features")
    print(f"  Train class distribution: 0={int((y_train==0).sum())}, 1={int((y_train==1).sum())}")

    return X_train, X_test, y_train, y_test, gene_names


def build_models(n_pos, n_neg):
    """
    Build the 4 ensemble classifiers.

    Per scikit-learn skill best practices:
    - Use class_weight='balanced' for imbalanced data
    - Set random_state for reproducibility
    - Use stratified CV
    """
    scale_pos_weight = n_neg / n_pos  # For XGBoost (ratio of neg/pos)

    models = {
        "RandomForest": RandomForestClassifier(
            n_estimators=200,
            max_depth=None,
            min_samples_split=5,
            min_samples_leaf=2,
            class_weight="balanced",
            random_state=RANDOM_STATE,
            n_jobs=N_JOBS,
        ),
        "XGBoost": XGBClassifier(
            n_estimators=200,
            max_depth=6,
            learning_rate=0.1,
            scale_pos_weight=scale_pos_weight,
            eval_metric="logloss",
            random_state=RANDOM_STATE,
            n_jobs=N_JOBS,
            verbosity=0,
        ),
        "LightGBM": LGBMClassifier(
            n_estimators=200,
            max_depth=-1,
            learning_rate=0.1,
            is_unbalance=True,
            random_state=RANDOM_STATE,
            n_jobs=N_JOBS,
            verbose=-1,
        ),
    }

    # Stacking meta-classifier: combines RF + XGB + LGBM predictions
    # using LogisticRegression as the final meta-learner
    # Per ml-engineer skill: stacking improves robustness
    #
    # IMPORTANT: Use separate model instances with n_jobs=1 for Stacking
    # to avoid nested parallelism (base models spawning worker pools
    # inside Stacking's sequential loop causes severe CPU contention)
    stacking_estimators = [
        ("rf", RandomForestClassifier(
            n_estimators=200, max_depth=None, min_samples_split=5,
            min_samples_leaf=2, class_weight="balanced",
            random_state=RANDOM_STATE, n_jobs=1,
        )),
        ("xgb", XGBClassifier(
            n_estimators=200, max_depth=6, learning_rate=0.1,
            scale_pos_weight=scale_pos_weight, eval_metric="logloss",
            random_state=RANDOM_STATE, n_jobs=1, verbosity=0,
        )),
        ("lgbm", LGBMClassifier(
            n_estimators=200, max_depth=-1, learning_rate=0.1,
            is_unbalance=True, random_state=RANDOM_STATE,
            n_jobs=1, verbose=-1,
        )),
    ]

    models["StackingEnsemble"] = StackingClassifier(
        estimators=stacking_estimators,
        final_estimator=LogisticRegression(
            class_weight="balanced",
            max_iter=1000,
            random_state=RANDOM_STATE,
        ),
        cv=StratifiedKFold(n_splits=3, shuffle=True, random_state=RANDOM_STATE),
        stack_method="predict_proba",
        n_jobs=1,
    )

    return models


def cross_validate_models(models, X_train, y_train):
    """Run stratified K-fold CV on each model and collect metrics."""
    print("\n" + "=" * 60)
    print(f"Cross-Validation ({N_CV_FOLDS}-Fold Stratified)")
    print("=" * 60)

    cv = StratifiedKFold(n_splits=N_CV_FOLDS, shuffle=True,
                         random_state=RANDOM_STATE)

    scoring = {
        "accuracy": "accuracy",
        "f1": "f1",
        "precision": "precision",
        "recall": "recall",
        "roc_auc": "roc_auc",
    }

    cv_results = {}

    for name, model in models.items():
        print(f"\n  Training {name}...")
        t0 = time.time()

        scores = cross_validate(
            model, X_train, y_train,
            cv=cv,
            scoring=scoring,
            return_train_score=False,
            n_jobs=1 if name == "StackingEnsemble" else N_JOBS,
        )

        elapsed = time.time() - t0

        result = {}
        for metric_name in scoring:
            key = f"test_{metric_name}"
            vals = scores[key]
            result[f"{metric_name}_mean"] = float(np.mean(vals))
            result[f"{metric_name}_std"] = float(np.std(vals))

        result["cv_time_seconds"] = elapsed
        cv_results[name] = result

        print(f"    -> AUC-ROC: {result['roc_auc_mean']:.4f} +/- {result['roc_auc_std']:.4f}")
        print(f"    -> F1:      {result['f1_mean']:.4f} +/- {result['f1_std']:.4f}")
        print(f"    -> Time:    {elapsed:.1f}s")

    return cv_results


def train_final_models(models, X_train, y_train):
    """Train each model on the full training set."""
    print("\n" + "=" * 60)
    print("Training Final Models on Full Training Set")
    print("=" * 60)

    trained_models = {}
    for name, model in models.items():
        print(f"  Training {name}...")
        t0 = time.time()
        model.fit(X_train, y_train)
        elapsed = time.time() - t0
        print(f"    -> Done in {elapsed:.1f}s")
        trained_models[name] = model

    return trained_models


def evaluate_on_test(trained_models, X_test, y_test):
    """Evaluate all trained models on the held-out test set."""
    print("\n" + "=" * 60)
    print("Test Set Evaluation")
    print("=" * 60)

    test_results = {}
    predictions = {}

    for name, model in trained_models.items():
        y_pred = model.predict(X_test)
        y_proba = model.predict_proba(X_test)[:, 1]

        metrics = {
            "accuracy": float(accuracy_score(y_test, y_pred)),
            "f1": float(f1_score(y_test, y_pred)),
            "precision": float(precision_score(y_test, y_pred)),
            "recall": float(recall_score(y_test, y_pred)),
            "roc_auc": float(roc_auc_score(y_test, y_proba)),
        }
        test_results[name] = metrics
        predictions[name] = {"y_pred": y_pred, "y_proba": y_proba}

        print(f"\n  {name}:")
        print(f"    Accuracy:  {metrics['accuracy']:.4f}")
        print(f"    F1 Score:  {metrics['f1']:.4f}")
        print(f"    Precision: {metrics['precision']:.4f}")
        print(f"    Recall:    {metrics['recall']:.4f}")
        print(f"    AUC-ROC:   {metrics['roc_auc']:.4f}")

    return test_results, predictions


def find_best_model(cv_results, test_results):
    """Select the best model based on test AUC-ROC."""
    best_name = max(test_results, key=lambda k: test_results[k]["roc_auc"])
    print(f"\n*** Best model: {best_name} "
          f"(Test AUC-ROC = {test_results[best_name]['roc_auc']:.4f})")
    return best_name


def plot_results(trained_models, predictions, test_results, y_test):
    """Generate publication-quality evaluation plots."""
    print("\n" + "=" * 60)
    print("Generating Evaluation Plots")
    print("=" * 60)

    os.makedirs(FIGURES_DIR, exist_ok=True)

    # --- 1. Model Comparison Bar Chart ---
    fig, ax = plt.subplots(figsize=(12, 6))
    metrics_to_plot = ["accuracy", "f1", "precision", "recall", "roc_auc"]
    model_names = list(test_results.keys())
    x = np.arange(len(model_names))
    width = 0.15

    colors = ["#2196F3", "#4CAF50", "#FF9800", "#9C27B0", "#F44336"]

    for i, metric in enumerate(metrics_to_plot):
        values = [test_results[m][metric] for m in model_names]
        bars = ax.bar(x + i * width, values, width, label=metric.upper(),
                      color=colors[i], edgecolor="white", linewidth=0.5)
        for bar, val in zip(bars, values):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.005,
                    f"{val:.3f}", ha="center", va="bottom", fontsize=7)

    ax.set_xlabel("Model", fontsize=12)
    ax.set_ylabel("Score", fontsize=12)
    ax.set_title("Ensemble ML Model Comparison -- TB vs. Control Classification",
                 fontsize=14, fontweight="bold")
    ax.set_xticks(x + width * 2)
    ax.set_xticklabels(model_names, fontsize=10)
    ax.legend(loc="lower right", fontsize=9)
    ax.set_ylim(0, 1.1)
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "model_comparison.png"), dpi=200)
    plt.close()

    # --- 2. ROC Curves ---
    fig, ax = plt.subplots(figsize=(8, 8))
    colors_roc = ["#1976D2", "#388E3C", "#F57C00", "#7B1FA2"]

    for (name, preds), color in zip(predictions.items(), colors_roc):
        RocCurveDisplay.from_predictions(
            y_test, preds["y_proba"],
            name=f"{name} (AUC={test_results[name]['roc_auc']:.3f})",
            ax=ax, color=color, linewidth=2,
        )

    ax.plot([0, 1], [0, 1], "k--", linewidth=1, alpha=0.5, label="Random (AUC=0.500)")
    ax.set_title("ROC Curves -- TB Classification Ensemble", fontsize=14, fontweight="bold")
    ax.legend(loc="lower right", fontsize=10)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "roc_curves.png"), dpi=200)
    plt.close()

    # --- 3. Confusion Matrices ---
    fig, axes = plt.subplots(1, 4, figsize=(20, 5))
    labels = ["Control", "TB-affected"]

    for ax, (name, preds) in zip(axes, predictions.items()):
        cm = confusion_matrix(y_test, preds["y_pred"])
        im = ax.imshow(cm, interpolation="nearest", cmap="Blues")
        ax.set_title(name, fontsize=12, fontweight="bold")

        # Add text annotations
        for i in range(2):
            for j in range(2):
                ax.text(j, i, f"{cm[i, j]:,}", ha="center", va="center",
                        fontsize=14, color="white" if cm[i, j] > cm.max() / 2 else "black")

        ax.set_xticks([0, 1])
        ax.set_yticks([0, 1])
        ax.set_xticklabels(labels, fontsize=9)
        ax.set_yticklabels(labels, fontsize=9)
        ax.set_ylabel("True Label")
        ax.set_xlabel("Predicted Label")

    fig.suptitle("Confusion Matrices -- Test Set", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "confusion_matrices.png"), dpi=200)
    plt.close()

    print(f"  -> Plots saved to '{FIGURES_DIR}/'")


def save_models_and_results(trained_models, cv_results, test_results, best_model_name):
    """Save trained models and results."""
    print("\n" + "=" * 60)
    print("Saving Models & Results")
    print("=" * 60)

    os.makedirs(MODELS_DIR, exist_ok=True)
    os.makedirs(RESULTS_DIR, exist_ok=True)

    # Save each trained model
    for name, model in trained_models.items():
        path = os.path.join(MODELS_DIR, f"{name}.joblib")
        joblib.dump(model, path)
        print(f"  -> Saved {name} --> {path}")

    # Save the best model name
    with open(os.path.join(MODELS_DIR, "best_model.txt"), "w") as f:
        f.write(best_model_name)

    # Save CV results
    cv_df = pd.DataFrame(cv_results).T
    cv_df.to_csv(os.path.join(RESULTS_DIR, "cv_results.csv"))
    print(f"  -> CV results --> {RESULTS_DIR}/cv_results.csv")

    # Save test results
    test_df = pd.DataFrame(test_results).T
    test_df.to_csv(os.path.join(RESULTS_DIR, "test_results.csv"))
    print(f"  -> Test results --> {RESULTS_DIR}/test_results.csv")

    # Combined summary
    summary = {
        "best_model": best_model_name,
        "best_test_auc": test_results[best_model_name]["roc_auc"],
        "best_test_f1": test_results[best_model_name]["f1"],
        "n_models": len(trained_models),
        "cv_folds": N_CV_FOLDS,
    }
    with open(os.path.join(RESULTS_DIR, "ml_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\n[OK] All models and results saved!")


if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    # Load data
    X_train, X_test, y_train, y_test, gene_names = load_data()

    # Calculate class balance for model configuration
    n_pos = int((y_train == 1).sum())
    n_neg = int((y_train == 0).sum())

    # Build models
    models = build_models(n_pos, n_neg)

    # Cross-validate
    cv_results = cross_validate_models(models, X_train, y_train)

    # Re-build models (CV mutates some) and train on full training set
    models = build_models(n_pos, n_neg)
    trained_models = train_final_models(models, X_train, y_train)

    # Evaluate on test set
    test_results, predictions = evaluate_on_test(trained_models, X_test, y_test)

    # Find best model
    best_model_name = find_best_model(cv_results, test_results)

    # Plot results
    plot_results(trained_models, predictions, test_results, y_test)

    # Save everything
    save_models_and_results(trained_models, cv_results, test_results, best_model_name)

    # Print final classification reports
    print("\n" + "=" * 60)
    print("Detailed Classification Reports (Test Set)")
    print("=" * 60)
    for name, preds in predictions.items():
        print(f"\n--- {name} ---")
        print(classification_report(y_test, preds["y_pred"],
                                    target_names=["Control", "TB-affected"]))

    print("\n--> Next: Run 04_model_improvement.py for hyperparameter tuning,")
    print("          or skip to 05_shap_analysis.py for disease signature extraction.")
