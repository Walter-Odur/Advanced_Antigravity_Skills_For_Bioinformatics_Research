# REINVENT4 Scoring Components Reference

Quick reference for all scoring components available in REINVENT4.
See `REINVENT4/configs/SCORING.md` for the full upstream documentation.

## Basic Molecular Properties (RDKit, no parameters)

| Component | Description | Recommended Transform |
|---|---|---|
| `Qed` | Drug-likeness (0–1) | None (already normalized) |
| `MolecularWeight` | Molecular weight | `double_sigmoid` low=200, high=500 |
| `TPSA` | Topological polar surface area | `double_sigmoid` low=20, high=140 |
| `SlogP` | LogP (Crippen) | `reverse_sigmoid` high=5, low=1 |
| `HBondAcceptors` | H-bond acceptors | `reverse_sigmoid` high=10, low=8 |
| `HBondDonors` | H-bond donors | `reverse_sigmoid` high=5, low=3 |
| `NumRotBond` | Rotatable bonds | `reverse_sigmoid` high=10, low=3 |
| `Csp3` | Fraction sp3 carbons | `sigmoid` low=0.2, high=0.5 |
| `NumHeavyAtoms` | Heavy atom count | `double_sigmoid` low=15, high=40 |
| `NumRings` | Ring count | `step` low=1, high=5 |
| `NumAromaticRings` | Aromatic ring count | `step` low=1, high=3 |
| `GraphLength` | Longest path in graph | `reverse_sigmoid` high=50, low=20 |
| `LargestRingSize` | Largest ring size | `step` low=5, high=7 |

## Drug-likeness & Synthesizability

| Component | Description | Parameters | Recommended Transform |
|---|---|---|---|
| `SAScore` | Synthetic accessibility (1–10, lower=easier) | None | `reverse_sigmoid` high=5, low=1 |
| `CustomAlerts` | SMARTS filter (0=match, 1=clean) | `smarts` (list) | None (binary filter) |

## Similarity & Cheminformatics

| Component | Description | Key Parameters |
|---|---|---|
| `TanimotoSimilarity` | Morgan FP similarity | `smiles`, `radius`, `use_counts` |
| `MatchingSubstructure` | Substructure match penalty | `smarts`, `use_chirality` |
| `GroupCount` | SMARTS pattern count | `smarts` |
| `MMP` | Matched molecular pairs | `reference_smiles`, `num_of_cuts` |

## External Scoring

| Component | Description | When to Use |
|---|---|---|
| `DockStream` | Docking (Vina, Glide, etc.) | Structure-based optimization |
| `Maize` | Workflow manager | Complex multi-step scoring |
| `ExternalProcess` | Custom executable | Any external scorer |
| `REST` | REST API scorer | Remote scoring servers |
| `ChemProp` / `ChemProp2` | D-MPNN QSAR models | Activity prediction |

## Transforms

| Transform | Shape | Key Params | Use When |
|---|---|---|---|
| `sigmoid` | S-curve up | `low`, `high`, `k` | Prefer values above threshold |
| `reverse_sigmoid` | S-curve down | `low`, `high`, `k` | Penalize values above threshold |
| `double_sigmoid` | Bell curve | `low`, `high`, `coef_div`, `coef_si`, `coef_se` | Target a specific range |
| `step` | Binary in-range | `low`, `high` | Hard cutoff |
| `right_step` | Binary ≥ threshold | `high` | Hard minimum |
| `left_step` | Binary ≤ threshold | `low` | Hard maximum |
| `exponential_decay` | Decay curve | `k` | Penalize increasing values |

## Aggregation Functions

| Function | TOML Key | Description |
|---|---|---|
| Geometric mean | `geometric_mean` | Default. Penalizes components that score zero. |
| Arithmetic mean | `arithmetic_mean` | More tolerant of individual low scores. |

## Common Scoring Recipes

### Lipinski-compliant drug

```toml
[[stage.scoring.component]]
[stage.scoring.component.Qed]
[[stage.scoring.component.Qed.endpoint]]
name = "QED"
weight = 1.0

[[stage.scoring.component]]
[stage.scoring.component.MolecularWeight]
[[stage.scoring.component.MolecularWeight.endpoint]]
name = "MW"
weight = 0.5
transform.type = "double_sigmoid"
transform.high = 500.0
transform.low = 200.0
transform.coef_div = 500.0
transform.coef_si = 20.0
transform.coef_se = 20.0
```

### Similarity to a known drug + novelty

```toml
[[stage.scoring.component]]
[stage.scoring.component.TanimotoSimilarity]
[[stage.scoring.component.TanimotoSimilarity.endpoint]]
name = "Similarity to erlotinib"
weight = 0.7
params.smiles = ["COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1"]
params.radius = [3]
params.use_counts = [true]
transform.type = "double_sigmoid"
transform.high = 0.6
transform.low = 0.3
transform.coef_div = 100.0
transform.coef_si = 20.0
transform.coef_se = 20.0
```
