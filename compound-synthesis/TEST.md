# TEST.md — Compound Synthesis Skill Test Suite

> Test suite for the `compound-synthesis` Antigravity skill (REINVENT4 wrapper).
> **All tests use minimal parameters to run in under 2 minutes total.**
> Tests marked ⚡ run in <5 seconds. Tests marked 🔬 run in <30 seconds.

---

## Setup

```bash
# All commands run from: e:\ANTIGRAVITY_WORKSHOP
# CLI script: .agents\skills\compound-synthesis\scripts\compound_synthesis.py
# Python env: .venv\Scripts\python.exe
# Output dir: compound_synthesis_tests\
```

---

## Test 1 ⚡ — Check Setup (validates installation)

```bash
python compound_synthesis.py check-setup \
    --reinvent-dir E:\ANTIGRAVITY_WORKSHOP\REINVENT4 \
    --output compound_synthesis_tests/t01_setup.json
```

**Expected output**:
- `python_version`: 3.11+
- `rdkit_available`: true
- `torch_available`: true or false (depends on env)
- `reinvent_importable`: true (if torch installed) or false
- `priors_found`: list of available .prior files
- `priors_missing`: list of missing .prior files
- Exit code: 0

**Pass criteria**: JSON written, no crash. Reports correct status for each dependency.

---

## Test 2 ⚡ — Generate Config: Staged Learning (drug-like)

```bash
python compound_synthesis.py generate-config \
    --mode staged_learning \
    --generator reinvent \
    --prior .reinvent \
    --scoring-profile drug-like \
    --num-steps 5 \
    --batch-size 16 \
    --device cpu \
    --output compound_synthesis_tests/t02_rl_config.toml
```

**Expected output**:
- Valid TOML file with `run_type = "staged_learning"`
- `device = "cpu"`
- `batch_size = 16`
- `max_steps = 5`
- 6 scoring components (Qed, MolecularWeight, SlogP, TPSA, NumRotBond, custom_alerts)
- `prior_file` path uses forward slashes (no backslashes)
- Exit code: 0

**Validation**:
```bash
python -c "import tomllib; d=tomllib.load(open('compound_synthesis_tests/t02_rl_config.toml','rb')); assert d['run_type']=='staged_learning'; assert d['parameters']['batch_size']==16; assert len(d['stage'][0]['scoring']['component'])==6; print('PASS')"
```

---

## Test 3 ⚡ — Generate Config: Transfer Learning

```bash
python compound_synthesis.py generate-config \
    --mode transfer_learning \
    --generator reinvent \
    --prior .reinvent \
    --smiles-file compound_synthesis_tests/t05_seeds.smi \
    --num-epochs 3 \
    --batch-size 16 \
    --device cpu \
    --output compound_synthesis_tests/t03_tl_config.toml
```

**Expected output**:
- Valid TOML with `run_type = "transfer_learning"`
- `num_epochs = 3`
- `batch_size = 16` (not 64 default)
- `input_model_file` points to reinvent.prior with forward slashes
- Exit code: 0

**Validation**:
```bash
python -c "import tomllib; d=tomllib.load(open('compound_synthesis_tests/t03_tl_config.toml','rb')); assert d['run_type']=='transfer_learning'; assert d['parameters']['num_epochs']==3; print('PASS')"
```

---

## Test 4 ⚡ — Generate Config: Sampling

```bash
python compound_synthesis.py generate-config \
    --mode sampling \
    --prior .reinvent \
    --num-smiles 50 \
    --device cpu \
    --output compound_synthesis_tests/t04_sampling_config.toml
```

**Expected output**:
- Valid TOML with `run_type = "sampling"`
- `num_smiles = 50`
- `model_file` path with forward slashes
- Exit code: 0

---

## Test 5 ⚡ — Prepare Seeds (validation + dedup)

```bash
python compound_synthesis.py prepare-seeds \
    --smiles "CCO" "c1ccccc1" "CC(=O)Oc1ccccc1C(=O)O" "INVALID_SMILES" "CCO" \
    --report compound_synthesis_tests/t05_seeds_report.json \
    --output compound_synthesis_tests/t05_seeds.smi
```

**Expected output**:
- `valid`: 3 (ethanol, benzene, aspirin)
- `invalid`: 1 ("INVALID_SMILES")
- `duplicates_removed`: 1 (second "CCO")
- Output `.smi` file has exactly 3 lines with canonical SMILES
- Exit code: 0

**Validation**:
```bash
python -c "lines=open('compound_synthesis_tests/t05_seeds.smi').readlines(); assert len(lines)==3, f'Expected 3, got {len(lines)}'; print('PASS')"
```

---

## Test 6 ⚡ — Generate Config: Scoring Mode

```bash
python compound_synthesis.py generate-config \
    --mode scoring \
    --smiles-file compound_synthesis_tests/t05_seeds.smi \
    --scoring-profile lead-like \
    --output compound_synthesis_tests/t06_scoring_config.toml
```

**Expected output**:
- Valid TOML with `run_type = "scoring"`
- `smiles_file` points to the seeds file
- Scoring section has 7 components (Qed, MW, SlogP, HBD, HBA, RotBonds, custom_alerts)
- Exit code: 0

---

## Test 7 ⚡ — Generate Config: Custom Components

```bash
python compound_synthesis.py generate-config \
    --mode staged_learning \
    --generator reinvent \
    --prior .reinvent \
    --component "Qed:weight=1.0" \
    --component "MolecularWeight:weight=0.5,transform=double_sigmoid,low=300,high=500" \
    --component "SAScore:weight=0.4,transform=reverse_sigmoid,high=4,low=1,k=0.5" \
    --num-steps 5 \
    --batch-size 16 \
    --device cpu \
    --output compound_synthesis_tests/t07_custom_config.toml
```

**Expected output**:
- Valid TOML with exactly 3 scoring components
- Component names: Qed, MolecularWeight, SAScore
- MolecularWeight has `transform.type = "double_sigmoid"` with `low = 300.0`
- Exit code: 0

---

## Test 8 ⚡ — Generate Config: Anti-TB Profile with Seeds

```bash
python compound_synthesis.py generate-config \
    --mode staged_learning \
    --generator reinvent \
    --prior .reinvent \
    --scoring-profile anti-tb \
    --seeds "Cc1ccnc(NC(=O)c2ccsc2NC(=O)c2ccn3cnnc3c2)n1" \
    --num-steps 5 \
    --batch-size 16 \
    --device cpu \
    --output compound_synthesis_tests/t08_antitb_config.toml
```

**Expected output**:
- Valid TOML with 5 scoring components (anti-tb profile)
- `smiles_file` auto-created at `t08_antitb_config_seeds.smi`
- Seed file contains 1 SMILES line
- Exit code: 0

---

## Test 9 ⚡ — Generate HPC Script

```bash
python compound_synthesis.py generate-hpc-script \
    --config compound_synthesis_tests/t02_rl_config.toml \
    --gpu a100 --ngpu 1 --mem 64G --time "24:00:00" \
    --conda-env reinvent4 \
    --job-name hydra_test \
    --output compound_synthesis_tests/t09_hpc.sh
```

**Expected output**:
- Valid bash script with SLURM directives
- Contains: `#SBATCH --gres=gpu:a100:1`
- Contains: `#SBATCH --mem=64G`
- Contains: `#SBATCH --time=24:00:00`
- Contains: `conda activate reinvent4`
- Contains: `reinvent` CLI command referencing the config
- Exit code: 0

**Validation**:
```bash
python -c "s=open('compound_synthesis_tests/t09_hpc.sh').read(); assert '--gres=gpu:a100:1' in s; assert 'conda activate reinvent4' in s; assert 'hydra_test' in s; print('PASS')"
```

---

## Test 10 ⚡ — Smart Device Detection (auto-detect)

```bash
python compound_synthesis.py generate-config \
    --mode staged_learning \
    --generator reinvent \
    --prior .reinvent \
    --scoring-profile drug-like \
    --num-steps 300 \
    --batch-size 128 \
    --output compound_synthesis_tests/t10_smart.toml
```

**Expected output** (no `--device` flag):
- Prints "No GPU detected, using CPU" (if no GPU) or "GPU detected: ..." (if GPU)
- Prints "Auto-selected device: cpu" or "cuda:0"
- If CPU + 300 steps: prints `[HPC RECOMMENDED]` and auto-generates `t10_smart_hpc.sh`
- Valid TOML with auto-detected device
- Exit code: 0

---

## Test 11 ⚡ — Analyze Results (mock CSV)

First, create a mock REINVENT output CSV:

```bash
python -c "
import csv, os
os.makedirs('compound_synthesis_tests', exist_ok=True)
with open('compound_synthesis_tests/mock_results.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['step', 'SMILES', 'total_score', 'Qed', 'SlogP'])
    w.writerow([1, 'CCO', 0.45, 0.32, 0.8])
    w.writerow([1, 'c1ccccc1', 0.62, 0.44, 0.7])
    w.writerow([2, 'CC(=O)O', 0.71, 0.55, 0.9])
    w.writerow([2, 'c1ccc(O)cc1', 0.83, 0.68, 0.85])
    w.writerow([3, 'CC(C)CC', 0.39, 0.28, 0.6])
print('Mock CSV created')
"
```

Then analyze:

```bash
python compound_synthesis.py analyze-results \
    --csv compound_synthesis_tests/mock_results.csv \
    --top-n 3 \
    --sort-by total_score \
    --output compound_synthesis_tests/t11_analysis.json
```

**Expected output**:
- `total_compounds`: 5
- `unique_molecules`: 5
- `max_score`: 0.83
- `mean_score`: 0.6 (approx)
- Top 3 compounds sorted by total_score descending:
  1. `c1ccc(O)cc1` (0.83)
  2. `CC(=O)O` (0.71)
  3. `c1ccccc1` (0.62)
- Exit code: 0

---

## Test 12 ⚡ — Download Priors (structure test)

```bash
python compound_synthesis.py download-priors \
    --reinvent-dir E:\ANTIGRAVITY_WORKSHOP\REINVENT4 \
    --models reinvent \
    --output compound_synthesis_tests/t12_download.json
```

**Expected output**:
- Attempts to fetch Zenodo record
- Either "Downloaded" or "already exists" or "error" (network dependent)
- JSON report written with keys: `priors_dir`, `downloaded`, `skipped`, `errors`
- Exit code: 0

---

## Full Test Runner Script

Run all tests in sequence (<2 minutes total):

```bash
mkdir -p compound_synthesis_tests

# Create mock data first
python -c "
import csv, os
os.makedirs('compound_synthesis_tests', exist_ok=True)
with open('compound_synthesis_tests/mock_results.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['step', 'SMILES', 'total_score', 'Qed', 'SlogP'])
    w.writerow([1, 'CCO', 0.45, 0.32, 0.8])
    w.writerow([1, 'c1ccccc1', 0.62, 0.44, 0.7])
    w.writerow([2, 'CC(=O)O', 0.71, 0.55, 0.9])
    w.writerow([2, 'c1ccc(O)cc1', 0.83, 0.68, 0.85])
    w.writerow([3, 'CC(C)CC', 0.39, 0.28, 0.6])
"

CS=".agents/skills/compound-synthesis/scripts/compound_synthesis.py"

echo "=== T01: check-setup ==="
python $CS check-setup --reinvent-dir E:\ANTIGRAVITY_WORKSHOP\REINVENT4 --output compound_synthesis_tests/t01_setup.json

echo "=== T05: prepare-seeds (run before T03) ==="
python $CS prepare-seeds --smiles "CCO" "c1ccccc1" "CC(=O)Oc1ccccc1C(=O)O" "INVALID" "CCO" --report compound_synthesis_tests/t05_report.json --output compound_synthesis_tests/t05_seeds.smi

echo "=== T02: generate-config staged_learning ==="
python $CS generate-config --mode staged_learning --generator reinvent --prior .reinvent --scoring-profile drug-like --num-steps 5 --batch-size 16 --device cpu --output compound_synthesis_tests/t02_rl.toml

echo "=== T03: generate-config transfer_learning ==="
python $CS generate-config --mode transfer_learning --prior .reinvent --smiles-file compound_synthesis_tests/t05_seeds.smi --num-epochs 3 --batch-size 16 --device cpu --output compound_synthesis_tests/t03_tl.toml

echo "=== T04: generate-config sampling ==="
python $CS generate-config --mode sampling --prior .reinvent --num-smiles 50 --device cpu --output compound_synthesis_tests/t04_sampling.toml

echo "=== T06: generate-config scoring ==="
python $CS generate-config --mode scoring --smiles-file compound_synthesis_tests/t05_seeds.smi --scoring-profile lead-like --output compound_synthesis_tests/t06_scoring.toml

echo "=== T07: custom components ==="
python $CS generate-config --mode staged_learning --prior .reinvent --component "Qed:weight=1.0" --component "MolecularWeight:weight=0.5,transform=double_sigmoid,low=300,high=500" --num-steps 5 --batch-size 16 --device cpu --output compound_synthesis_tests/t07_custom.toml

echo "=== T08: anti-tb with seeds ==="
python $CS generate-config --mode staged_learning --prior .reinvent --scoring-profile anti-tb --seeds "Cc1ccnc(NC(=O)c2ccsc2NC(=O)c2ccn3cnnc3c2)n1" --num-steps 5 --batch-size 16 --device cpu --output compound_synthesis_tests/t08_antitb.toml

echo "=== T09: HPC script ==="
python $CS generate-hpc-script --config compound_synthesis_tests/t02_rl.toml --gpu a100 --ngpu 1 --mem 64G --time "24:00:00" --job-name hydra_test --output compound_synthesis_tests/t09_hpc.sh

echo "=== T10: smart device detection ==="
python $CS generate-config --mode staged_learning --prior .reinvent --scoring-profile drug-like --num-steps 300 --batch-size 128 --output compound_synthesis_tests/t10_smart.toml

echo "=== T11: analyze-results ==="
python $CS analyze-results --csv compound_synthesis_tests/mock_results.csv --top-n 3 --sort-by total_score --output compound_synthesis_tests/t11_analysis.json

echo "=== T12: download-priors ==="
python $CS download-priors --reinvent-dir E:\ANTIGRAVITY_WORKSHOP\REINVENT4 --models reinvent --output compound_synthesis_tests/t12_download.json

echo ""
echo "=== ALL TESTS COMPLETE ==="
```

---

## Expected Test Times

| Test | Subcommand | Expected Time |
|:----:|---|---|
| T01 | check-setup | <1s |
| T02 | generate-config (RL) | <1s |
| T03 | generate-config (TL) | <1s |
| T04 | generate-config (sampling) | <1s |
| T05 | prepare-seeds | <1s |
| T06 | generate-config (scoring) | <1s |
| T07 | generate-config (custom) | <1s |
| T08 | generate-config (anti-tb) | <1s |
| T09 | generate-hpc-script | <1s |
| T10 | smart device detection | <1s |
| T11 | analyze-results | <1s |
| T12 | download-priors | <30s (network) |
| **Total** | | **<45 seconds** |

---

## Notes

- Tests T01–T11 are **fully offline** (no network needed)
- Test T12 requires network access to Zenodo
- All tests use `--device cpu` and minimal parameters (5 steps, batch=16)
- No test requires PyTorch or REINVENT4 to be installed (except T12 download)
- TOML validation can be done with Python's built-in `tomllib` (3.11+)
- Tests generate ~15 output files in `compound_synthesis_tests/`
