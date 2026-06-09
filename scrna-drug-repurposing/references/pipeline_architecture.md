# Pipeline Architecture Reference

## Data Flow

```
Raw CSV files ──► AnnData (.h5ad) ──► NumPy arrays (.npz) ──► Trained models (.joblib)
                                                              ──► SHAP values
                                                              ──► Disease signature (.csv)
                                                              ──► Drug candidates (.csv)
                                                              ──► Figures (.png/.pdf)
```

## File I/O Contracts

### Step 1 → Step 2
| Output | Format | Description |
|---|---|---|
| `data/processed_tb_lung.h5ad` | AnnData HDF5 | QC'd, normalized, HVG-filtered. `.raw` stores all genes. |

### Step 2 → Step 3
| Output | Format | Description |
|---|---|---|
| `data/ml_ready/train_data.npz` | Compressed NumPy | Keys: `X` (float32 matrix), `y` (int labels) |
| `data/ml_ready/test_data.npz` | Compressed NumPy | Same structure as train |
| `data/ml_ready/gene_names.csv` | CSV | Single column `gene` — feature names for SHAP |
| `data/ml_ready/train_cells.csv` | CSV | Cell barcode IDs for traceability |
| `data/ml_ready/test_cells.csv` | CSV | Cell barcode IDs for traceability |
| `data/ml_ready/data_summary.json` | JSON | n_train, n_test, n_features, class counts |

### Step 3 → Steps 4, 5
| Output | Format | Description |
|---|---|---|
| `models/<ModelName>.joblib` | Joblib | Serialized trained scikit-learn model |
| `models/best_model.txt` | Text | Name of the best model (e.g., "LightGBM") |
| `results/cv_results.csv` | CSV | Cross-validation metrics per model |
| `results/test_results.csv` | CSV | Test set metrics per model |
| `results/ml_summary.json` | JSON | Best model name, AUC, F1, config |

### Step 5 → Step 6
| Output | Format | Description |
|---|---|---|
| `results/disease_signature.csv` | CSV | Top N genes: gene, mean_abs_shap, mean_shap, direction, rank |
| `results/full_gene_ranking.csv` | CSV | All genes ranked by SHAP importance |
| `results/signature_UP_genes.txt` | Text | One gene per line — upregulated in disease |
| `results/signature_DOWN_genes.txt` | Text | One gene per line — downregulated in disease |
| `results/shap_summary.json` | JSON | n_signature_genes, n_up, n_down, top_10_genes |

### Step 6 → Figures
| Output | Format | Description |
|---|---|---|
| `results/repurposed_drugs.csv` | CSV | Full drug repurposing results with scores |
| `results/drug_repurposing_summary.json` | JSON | Summary statistics |

## Directory Structure Convention

```
project_root/
├── 01_load_and_qc.py          # Numbered for execution order
├── 02_prepare_ml_data.py
├── 03_ensemble_ml.py
├── 04_model_improvement.py
├── 05_shap_analysis.py
├── 06_drug_repurposing.py
├── generate_manuscript_figures.py
├── data/                       # Input data + processed intermediates
│   ├── raw_counts.csv
│   ├── metadata.csv
│   ├── processed.h5ad          # Step 1 output (cached)
│   └── ml_ready/               # Step 2 output (cached)
├── models/                     # Step 3 output (cached)
├── results/                    # Steps 3-6 outputs
│   └── figures/
│       ├── qc/
│       ├── ml/
│       ├── shap/
│       ├── tuning/
│       ├── repurposing/
│       └── manuscript/         # Publication figures
└── manuscript/                 # LaTeX/Word manuscript
```

## Compute & Memory Requirements

### Estimated Resources by Dataset Size

| Dataset Size | Step 1 (QC) | Step 2 (Split) | Step 3 (ML Train) | Step 5 (SHAP) | Step 6 (Drugs) | Total |
|---|---|---|---|---|---|---|
| **5K cells × 2K genes** | ~30s, 1 GB | ~5s, 0.5 GB | ~2 min, 1 GB | ~1 min, 1 GB | ~30s, 0.5 GB | **~5 min** |
| **50K cells × 3K genes** | ~2 min, 4 GB | ~15s, 2 GB | ~15 min, 4 GB | ~10 min, 3 GB | ~1 min, 0.5 GB | **~30 min** |
| **200K cells × 5K genes** | ~10 min, 16 GB | ~1 min, 8 GB | ~1 hr, 12 GB | ~30 min, 8 GB | ~2 min, 0.5 GB | **~2 hrs** |
| **1M cells × 5K genes** | ~1 hr, 32 GB | ~5 min, 16 GB | ~4 hrs, 24 GB | ~2 hrs, 16 GB | ~5 min, 0.5 GB | **~8 hrs** |

> Memory = peak RSS. Times on 8-core CPU. SHAP uses subsampling (2,000 cells) for datasets > 2,000 cells.

### GPU Requirements (Optional Advanced Steps)

| Tool | GPU Needed? | Purpose |
|---|---|---|
| CellBender | ✅ Recommended (CUDA) | Ambient RNA removal (deep generative model) |
| scVI | ✅ Recommended (CUDA) | Batch correction / integration |
| SOLO | ✅ Recommended (CUDA) | Deep learning doublet detection |
| TreeExplainer | ❌ CPU only | SHAP computation |
| Ensemble ML | ❌ CPU only | RF/XGBoost/LightGBM training |

### Parallelization Strategy

| Step | Parallelism | Notes |
|---|---|---|
| ML training (standalone) | `n_jobs=-1` | Use all CPU cores |
| ML training (inside Stacking) | `n_jobs=1` | Avoid nested parallelism deadlock on Windows |
| Cross-validation | `n_jobs=-1` for standalone; `n_jobs=1` for Stacking | Same deadlock avoidance |
| SHAP TreeExplainer | Single-threaded | TreeExplainer is inherently sequential |
| Scanpy preprocessing | Automatic | Uses internal parallelism |

### Large Dataset Strategies (>100K cells)

- Use `sc.read_h5ad(filename, backed='r')` for memory-mapped reading
- Subsample balanced training sets when full dataset exceeds 200K cells
- Use AnnData `.chunked_X()` for out-of-memory iteration
- Consider Dask arrays for datasets exceeding available RAM

## Step Dependency Graph

```
01_load_and_qc.py ──► 02_prepare_ml_data.py ──► 03_ensemble_ml.py ──┬──► 04_model_improvement.py ──┐
                                                                     │                              │
                                                                     └──► 05_shap_analysis.py ◄─────┘
                                                                                   │
                                                                                   ▼
                                                                          06_drug_repurposing.py
                                                                                   │
                                                                                   ▼
                                                                     generate_manuscript_figures.py
                                                                     (reads all results/* outputs)
```

**Note:** Step 04 (model improvement/tuning) is optional. If skipped, Step 05 uses the model from Step 03 directly.

## Caching Strategy

Each script checks for its output files before running:

```python
if os.path.exists(OUTPUT_FILE) and not FORCE_RERUN:
    print(f"Found cached {OUTPUT_FILE}, skipping...")
    return load_cached(OUTPUT_FILE)
```

Set `FORCE_RERUN = True` at the top of any script to regenerate its outputs.

## Reproducibility Checklist

- [ ] `RANDOM_STATE = 42` used in all scripts
- [ ] `train_test_split(..., stratify=y, random_state=42)` for deterministic splits
- [ ] `StratifiedKFold(..., shuffle=True, random_state=42)` for CV
- [ ] NumPy/SHAP subsampling uses `np.random.RandomState(42)`
- [ ] All models set `random_state=42`
- [ ] Gene names saved alongside feature matrices
- [ ] Cell IDs saved for full traceability back to AnnData
