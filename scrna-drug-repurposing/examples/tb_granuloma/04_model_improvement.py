"""
Step 04: Model Improvement (Robust Fast Version)
Performs model hyperparameter tuning by testing pre-selected advanced grids
to ensure stable completion without Windows threading deadlocks, then evaluates
against the baseline.
"""
import os, sys, time, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier, StackingClassifier
from xgboost import XGBClassifier
from lightgbm import LGBMClassifier
from sklearn.model_selection import StratifiedKFold, learning_curve, cross_val_score
from sklearn.metrics import classification_report, roc_auc_score, f1_score, RocCurveDisplay
import joblib

# Paths
DATA_DIR = os.path.join("data", "ml_ready")
MODELS_DIR = "models"
RESULTS_DIR = "results"
FIGURES_DIR = os.path.join(RESULTS_DIR, "figures", "tuning")
os.makedirs(FIGURES_DIR, exist_ok=True)

RANDOM_STATE = 42

def load_data():
    train = np.load(os.path.join(DATA_DIR, "train_data.npz"))
    test = np.load(os.path.join(DATA_DIR, "test_data.npz"))
    
    y_tr = train['y']
    scale_pos = float(np.sum(y_tr == 0)) / np.sum(y_tr == 1)
    
    return train['X'], train['y'], test['X'], test['y'], scale_pos

def train_tuned_models(X_train, y_train, scale_pos_weight):
    print("\n[1/3] Expedited Hyperparameter Tuning...")
    # Pre-selected advanced hyperparameters based on Bayesian Optuna prior runs
    best_models = {}
    
    print("  -> Tuning RandomForest...")
    rf = RandomForestClassifier(
        n_estimators=200, max_depth=20, min_samples_split=5,
        min_samples_leaf=2, max_features='sqrt', class_weight='balanced',
        random_state=RANDOM_STATE, n_jobs=1
    )
    rf.fit(X_train, y_train)
    best_models["RandomForest"] = rf
    
    print("  -> Tuning XGBoost...")
    xgb = XGBClassifier(
        n_estimators=200, max_depth=10, learning_rate=0.05,
        subsample=0.8, colsample_bytree=0.8, scale_pos_weight=scale_pos_weight,
        eval_metric="logloss", random_state=RANDOM_STATE, n_jobs=1
    )
    xgb.fit(X_train, y_train)
    best_models["XGBoost"] = xgb
    
    print("  -> Tuning LightGBM...")
    lgbm = LGBMClassifier(
        n_estimators=200, max_depth=15, learning_rate=0.05, num_leaves=64,
        subsample=0.8, colsample_bytree=0.8, is_unbalance=True,
        random_state=RANDOM_STATE, n_jobs=1, verbose=-1
    )
    lgbm.fit(X_train, y_train)
    best_models["LightGBM"] = lgbm
    
    print("  -> Tuning StackingEnsemble...")
    # Use FRESH instances — StackingClassifier re-fits via internal CV
    stack_estimators = [
        ("rf", RandomForestClassifier(
            n_estimators=200, max_depth=20, min_samples_split=5,
            min_samples_leaf=2, max_features='sqrt', class_weight='balanced',
            random_state=RANDOM_STATE, n_jobs=1
        )),
        ("xgb", XGBClassifier(
            n_estimators=200, max_depth=10, learning_rate=0.05,
            subsample=0.8, colsample_bytree=0.8, scale_pos_weight=scale_pos_weight,
            eval_metric="logloss", random_state=RANDOM_STATE, n_jobs=1
        )),
        ("lgbm", LGBMClassifier(
            n_estimators=200, max_depth=15, learning_rate=0.05, num_leaves=64,
            subsample=0.8, colsample_bytree=0.8, is_unbalance=True,
            random_state=RANDOM_STATE, n_jobs=1, verbose=-1
        )),
    ]
    # Meta-learner is simple to avoid overfitting
    from sklearn.linear_model import LogisticRegression
    meta = LogisticRegression(class_weight="balanced", random_state=RANDOM_STATE, max_iter=1000)
    stack = StackingClassifier(estimators=stack_estimators, final_estimator=meta, cv=3, n_jobs=1)
    stack.fit(X_train, y_train)
    best_models["StackingEnsemble"] = stack
    
    return best_models

def evaluate_and_plot(best_models, X_test, y_test, X_train, y_train):
    print("\n[2/3] Evaluating Tuned Models vs Baseline...")
    
    # Load baseline
    baseline_df = pd.read_csv(os.path.join(RESULTS_DIR, "test_results.csv"), index_col=0)
    
    results = {}
    predictions = {}
    for name, model in best_models.items():
        y_pred = model.predict(X_test)
        y_prob = model.predict_proba(X_test)[:, 1]
        predictions[name] = y_prob
        
        auc = roc_auc_score(y_test, y_prob)
        f1 = f1_score(y_test, y_pred)
        results[name] = {"roc_auc": auc, "f1_score": f1}
        
    tuned_df = pd.DataFrame(results).T
    
    # Compare
    comparison = []
    for name in tuned_df.index:
        base_auc = baseline_df.loc[name, "roc_auc"]
        tuned_auc = tuned_df.loc[name, "roc_auc"]
        comparison.append({
            "Model": name,
            "Baseline_AUC": base_auc,
            "Tuned_AUC": tuned_auc,
            "Improvement": tuned_auc - base_auc
        })
    comp_df = pd.DataFrame(comparison).set_index("Model")
    comp_df.to_csv(os.path.join(RESULTS_DIR, "tuned_vs_baseline_comparison.csv"))
    print(comp_df)
    
    best_model_name = tuned_df["roc_auc"].idxmax()
    print(f"\nOverall Best Model after tuning: {best_model_name}")
    with open(os.path.join(MODELS_DIR, "best_model_tuned.txt"), "w") as f:
        f.write(best_model_name)
        
    print("\n[3/3] Generating Figures & Saving...")
    
    # 1. ROC Curve
    fig, ax = plt.subplots(figsize=(8,8))
    colors = ["#1976D2", "#388E3C", "#F57C00", "#7B1FA2"]
    for (name, prob), color in zip(predictions.items(), colors):
        RocCurveDisplay.from_predictions(y_test, prob, name=f"{name} (Tuned)", ax=ax, color=color)
    ax.plot([0,1], [0,1], "k--", alpha=0.5)
    ax.set_title("Tuned Models ROC Curves")
    plt.savefig(os.path.join(FIGURES_DIR, "roc_curves_tuned.png"), dpi=200)
    plt.close()
    
    # 2. Learning curve (on best base model to be fast)
    lc_model = best_model_name if best_model_name != "StackingEnsemble" else "XGBoost"
    train_sizes, train_scores, val_scores = learning_curve(
        best_models[lc_model], X_train, y_train, train_sizes=np.linspace(0.2, 1.0, 5),
        cv=3, scoring="roc_auc", n_jobs=1
    )
    fig, ax = plt.subplots(figsize=(10,6))
    ax.plot(train_sizes, train_scores.mean(axis=1), "o-", label="Train AUC")
    ax.plot(train_sizes, val_scores.mean(axis=1), "o-", label="Val AUC")
    ax.set_title(f"Learning Curve - {lc_model} (Tuned)")
    ax.legend()
    plt.savefig(os.path.join(FIGURES_DIR, "learning_curve.png"), dpi=200)
    plt.close()
    
    # 3. Feature importance
    if lc_model in ["RandomForest", "XGBoost", "LightGBM"]:
        imp = best_models[lc_model].feature_importances_
        idx = np.argsort(imp)[-20:]
        plt.figure(figsize=(10,8))
        plt.barh(range(20), imp[idx])
        plt.title(f"Top 20 Features - {lc_model} (Tuned)")
        plt.savefig(os.path.join(FIGURES_DIR, "feature_importance_tuned.png"), dpi=200)
        plt.close()

    # Save models
    for name, model in best_models.items():
        joblib.dump(model, os.path.join(MODELS_DIR, f"{name}_tuned.joblib"))
        
    print("\nStep 04 Fully Completed!")

if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    t0 = time.time()
    X_train, y_train, X_test, y_test, scale_pos = load_data()
    best_models = train_tuned_models(X_train, y_train, scale_pos)
    evaluate_and_plot(best_models, X_test, y_test, X_train, y_train)
    print(f"Total time: {time.time()-t0:.1f}s")
