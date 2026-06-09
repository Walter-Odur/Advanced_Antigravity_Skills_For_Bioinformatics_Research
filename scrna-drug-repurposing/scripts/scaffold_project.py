"""
Scaffold a New scRNA-seq Drug Repurposing Project
====================================================
Creates the directory structure and starter scripts for a new project.

Usage:
    python scaffold_project.py --name my_disease_study --path /path/to/project
    python scaffold_project.py  # Creates in current directory with default name
"""

import os
import argparse
import shutil


def create_directory_structure(project_dir):
    """Create the standard project directory structure."""
    directories = [
        "data",
        "data/ml_ready",
        "models",
        "results",
        "results/figures",
        "results/figures/qc",
        "results/figures/ml",
        "results/figures/shap",
        "results/figures/tuning",
        "results/figures/repurposing",
        "results/figures/manuscript",
        "manuscript",
    ]

    for d in directories:
        path = os.path.join(project_dir, d)
        os.makedirs(path, exist_ok=True)
        print(f"  Created: {path}")


def create_gitignore(project_dir):
    """Create a .gitignore for the project."""
    gitignore_content = """# Data files (too large for git)
data/*.csv
data/*.h5ad
data/*.h5
data/ml_ready/*.npz

# Models
models/*.joblib

# Python
__pycache__/
*.pyc
.venv/
venv/

# IDE
.vscode/
.idea/

# OS
.DS_Store
Thumbs.db

# Results (regenerable)
results/figures/
"""
    with open(os.path.join(project_dir, ".gitignore"), "w") as f:
        f.write(gitignore_content)
    print(f"  Created: .gitignore")


def create_readme(project_dir, project_name):
    """Create a starter README."""
    readme = f"""# {project_name}

## Overview

Single-cell RNA-seq ensemble ML pipeline for disease signature extraction and
computational drug repurposing.

Built using the `scrna-drug-repurposing` Antigravity skill.

## Pipeline

```bash
python 01_load_and_qc.py
python 02_prepare_ml_data.py
python 03_ensemble_ml.py
python 04_model_improvement.py
python 05_shap_analysis.py
python 06_drug_repurposing.py
python generate_manuscript_figures.py
```

## Requirements

```bash
pip install scanpy scikit-learn xgboost lightgbm shap optuna matplotlib seaborn pandas numpy scipy requests gseapy joblib
```

## Data

Place your scRNA-seq data files in the `data/` directory:
- `raw_counts.csv` — Expression matrix
- `metadata.csv` — Cell annotations with disease status column

## Configuration

Edit the configuration block at the top of each script to match your dataset.
See `config_template.py` for all available parameters.
"""
    with open(os.path.join(project_dir, "README.md"), "w") as f:
        f.write(readme)
    print(f"  Created: README.md")


def create_placeholder_scripts(project_dir):
    """Create placeholder pipeline scripts with TODOs."""
    scripts = {
        "01_load_and_qc.py": '"""\nStep 01: Load Raw Data & Quality Control\n\nTODO: Adapt from the scrna-drug-repurposing skill examples.\nSee: examples/tb_granuloma/01_load_and_qc.py\n"""\n\nprint("Step 01: Not yet implemented. Use the skill reference to build this step.")\n',
        "02_prepare_ml_data.py": '"""\nStep 02: Prepare ML-Ready Data\n\nTODO: Adapt from the scrna-drug-repurposing skill examples.\nSee: examples/tb_granuloma/02_prepare_ml_data.py\n"""\n\nprint("Step 02: Not yet implemented.")\n',
        "03_ensemble_ml.py": '"""\nStep 03: Ensemble ML Training & Evaluation\n\nTODO: Adapt from the scrna-drug-repurposing skill examples.\nSee: examples/tb_granuloma/03_ensemble_ml.py\n"""\n\nprint("Step 03: Not yet implemented.")\n',
        "04_model_improvement.py": '"""\nStep 04: Model Improvement (Hyperparameter Tuning)\n\nTODO: Adapt from the scrna-drug-repurposing skill examples.\n"""\n\nprint("Step 04: Not yet implemented.")\n',
        "05_shap_analysis.py": '"""\nStep 05: SHAP Analysis & Disease Signature Extraction\n\nTODO: Adapt from the scrna-drug-repurposing skill examples.\n"""\n\nprint("Step 05: Not yet implemented.")\n',
        "06_drug_repurposing.py": '"""\nStep 06: Drug Repurposing via CMap/L1000\n\nTODO: Adapt from the scrna-drug-repurposing skill examples.\n"""\n\nprint("Step 06: Not yet implemented.")\n',
        "generate_manuscript_figures.py": '"""\nGenerate Publication-Quality Manuscript Figures\n\nTODO: Adapt from the scrna-drug-repurposing skill examples.\n"""\n\nprint("Figure generation: Not yet implemented.")\n',
    }

    for filename, content in scripts.items():
        filepath = os.path.join(project_dir, filename)
        with open(filepath, "w") as f:
            f.write(content)
        print(f"  Created: {filename}")


def main():
    parser = argparse.ArgumentParser(
        description="Scaffold a new scRNA-seq drug repurposing project"
    )
    parser.add_argument(
        "--name", default="scrna_drug_repurposing_project",
        help="Project name (used in README)"
    )
    parser.add_argument(
        "--path", default=".",
        help="Directory to create the project in (default: current directory)"
    )
    args = parser.parse_args()

    project_dir = os.path.abspath(args.path)

    print(f"\n{'=' * 60}")
    print(f"Scaffolding: {args.name}")
    print(f"Location:    {project_dir}")
    print(f"{'=' * 60}\n")

    create_directory_structure(project_dir)
    create_gitignore(project_dir)
    create_readme(project_dir, args.name)
    create_placeholder_scripts(project_dir)

    print(f"\n{'=' * 60}")
    print(f"Project scaffolded successfully!")
    print(f"\nNext steps:")
    print(f"  1. Place your data files in {os.path.join(project_dir, 'data')}/")
    print(f"  2. Edit the configuration in each script")
    print(f"  3. Run: python 01_load_and_qc.py")
    print(f"{'=' * 60}\n")


if __name__ == "__main__":
    main()
