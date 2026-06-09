"""
Project Configuration Template
================================
Adapt this template for new scRNA-seq drug repurposing projects.
Copy this file to your project root and modify the values below.

Usage:
    from config_template import *
"""

# ============================================================
# PROJECT METADATA
# ============================================================
PROJECT_NAME = "my_disease_study"
DATA_SOURCE = "SCP_or_GEO_accession_number"
ORGANISM = "Human"  # "Human" or "Mouse"

# ============================================================
# FILE PATHS
# ============================================================
DATA_DIR = "data"
RAW_COUNTS_FILE = "raw_counts.csv"        # Expression matrix (cells x genes or genes x cells)
METADATA_FILE = "metadata.csv"            # Cell annotations with disease status
OUTPUT_H5AD = "processed.h5ad"            # Processed AnnData output

ML_READY_DIR = "data/ml_ready"
MODELS_DIR = "models"
RESULTS_DIR = "results"
FIGURES_DIR = "results/figures"

# ============================================================
# QC THRESHOLDS (Step 01)
# ============================================================
MIN_GENES_PER_CELL = 200       # Min genes detected per cell (removes empty droplets)
MIN_CELLS_PER_GENE = 3         # Min cells expressing a gene (removes rare genes)

# Mitochondrial threshold — adjust per tissue type:
#   PBMCs/Blood:    5-10%
#   Solid tissues:  15-25% (lung, liver, brain have higher MT content)
#   Tumors:         10-20% (variable)
MAX_PCT_MITO = 20

# HVG selection parameters (Scanpy defaults, good for most datasets)
HVG_MIN_MEAN = 0.0125
HVG_MAX_MEAN = 3
HVG_MIN_DISP = 0.5

# ============================================================
# LABEL MAPPING (Step 02)
# ============================================================
# Map your disease status column to binary labels
# Key = value in metadata Disease_Status column
# Value = 0 (control) or 1 (disease)
DISEASE_STATUS_COLUMN = "Disease_Status"

LABEL_MAP = {
    "Disease_Group": 1,
    "Control_Group": 0,
    # Add more groups as needed:
    # "Disease_Subtype": 1,
    # "Healthy": 0,
}

# Groups to EXCLUDE entirely (confounders, ambiguous, etc.)
EXCLUDED_GROUPS = []

# ============================================================
# ML CONFIGURATION (Steps 03-04)
# ============================================================
TEST_SIZE = 0.20               # Fraction held out for testing
N_CV_FOLDS = 5                 # Cross-validation folds
RANDOM_STATE = 42              # Global random seed for reproducibility
N_JOBS = -1                    # CPU cores for training (-1 = all available)

# Model hyperparameters (defaults from the TB pipeline)
RF_N_ESTIMATORS = 200
XGB_N_ESTIMATORS = 200
XGB_MAX_DEPTH = 6
XGB_LEARNING_RATE = 0.1
LGBM_N_ESTIMATORS = 200
LGBM_LEARNING_RATE = 0.1

# ============================================================
# SHAP CONFIGURATION (Step 05)
# ============================================================
TOP_N_GENES = 100              # Number of genes in the disease signature
SHAP_SAMPLE_SIZE = 2000        # Subsample size for SHAP computation
MAX_DISPLAY = 30               # Max genes shown in SHAP plots
TOP_N_DEPENDENCE = 5           # Number of dependence plots to generate

# ============================================================
# DRUG REPURPOSING (Step 06)
# ============================================================
# Gene set libraries to query via GSEApy
GENE_SET_LIBRARIES = [
    "L1000_Ligand_Pert_up",
    "L1000_Ligand_Pert_down",
    "Drug_Perturbations_from_GEO_up",
    "Drug_Perturbations_from_GEO_down",
]

# Significance threshold for drug candidates
DRUG_PVALUE_THRESHOLD = 0.05

# Minimum combined prioritization score to include in final candidates
MIN_COMBINED_SCORE = 1.0

# ============================================================
# FIGURE GENERATION (Step 07)
# ============================================================
DISEASE_COLOR = '#E63946'      # Color for disease class in plots
CONTROL_COLOR = '#457B9D'      # Color for control class in plots
FIGURE_DPI = 300               # Resolution for saved figures
SAVE_PDF = True                # Also save PDF versions for manuscript

# ============================================================
# CACHING
# ============================================================
FORCE_RERUN = False            # Set True to regenerate all outputs
