# ADMET Thresholds Reference

## Lipinski Rule of Five

| Property | Threshold | Description |
|----------|-----------|-------------|
| MW | <= 500 Da | Molecular weight |
| LogP | <= 5 | Partition coefficient (lipophilicity) |
| HBD | <= 5 | Hydrogen bond donors (OH + NH count) |
| HBA | <= 10 | Hydrogen bond acceptors (N + O count) |

**Pass criteria**: <= 1 violation. Compounds with 0-1 violations are
considered drug-like.

**Note**: RDKit uses the Lipinski definition for HBA counting, which may differ
from some tools. HBA counts only N and O atoms (not S), and excludes certain
delocalized systems.

## Beyond Rule of Five (bRo5)

For larger molecules (e.g., macrocycles, PROTACs), relaxed thresholds apply:

| Property | Threshold |
|----------|-----------|
| MW | <= 1000 Da |
| LogP | <= 10 |
| HBD | <= 5 |
| HBA | <= 15 |
| PSA | <= 250 A^2 |

## PAINS Alerts

PAINS (Pan-Assay Interference Compounds) are frequent hitters in HTS. The
filter uses RDKit's FilterCatalog with three tiers:

| Catalog | # Filters | Description |
|---------|-----------|-------------|
| PAINS_A | 16 | Most problematic substructures |
| PAINS_B | 55 | Moderate concern |
| PAINS_C | 409 | Lower confidence alerts |

### Common PAINS Substructures

| Alert | SMARTS Pattern | Notes |
|-------|----------------|-------|
| Rhodanine | `[#6](=[#16])(-[#7])-[#16]` | rhod_sat_A(33) |
| Quinone | `O=C1C=CC(=O)C=C1` | Often reactive |
| Catechol | `c1cc(O)c(O)cc1` | Redox cycling |
| Hydroxyphenyl hydrazone | Multiple | Assay interference |

**Pass criteria**: Zero PAINS alerts across all three catalogs.

## Synthetic Accessibility (SA) Score

Ertl & Schuffenhauer score on a 1-10 scale:

| Range | Interpretation |
|-------|----------------|
| 1-3 | Easy to synthesize |
| 3-5 | Moderately difficult |
| 5-7 | Difficult |
| 7-10 | Very difficult / impractical |

**Default threshold**: <= 6.0 (configurable via `--sa-threshold`)

### Reference Values

| Compound | SA Score |
|----------|----------|
| Aspirin | ~1.6 |
| Ibuprofen | ~2.2 |
| Caffeine | ~2.3 |
| Taxol | ~7.9 |
| Vancomycin | ~8.5 |

## Cascade Order

The ADMET filter applies checks in this order:
1. **Lipinski** (fast, property-based)
2. **PAINS** (substructure matching)
3. **SA Score** (computational, slower)

A compound must pass ALL three to receive an overall PASS.

## Two-Tier Filtering Strategy

The wizard supports two strictness profiles selected via `--strictness`:

### Tier 1: Relaxed (post-docking → REINVENT seeds)

Applied after docking to keep more hits as training seeds for REINVENT.
Skips metabolism (CYP), most toxicity (hERG, DILI), and medicinal chem rules.

| Property | Relaxed Threshold | vs Strict |
|----------|-------------------|-----------|
| MW | ≤ 700 Da | +200 |
| LogP | ≤ 7.0 | +2.0 |
| HBD | ≤ 7 | +2 |
| HBA | ≤ 15 | +5 |
| TPSA | ≤ 200 Å² | +60 |
| LogS | ≥ -8.0 | -2.0 |
| PPB | ≤ 99% | +4% |
| PAINS | ≤ 1 alert | +1 |
| Brenk | ≤ 2 alerts | +2 |
| SA Score | ≤ 8.0 | +2.0 |
| CYP, hERG, DILI | **SKIPPED** | — |
| Lipinski/Pfizer rules | **SKIPPED** | — |

**14 properties checked** (vs 34 in strict mode).

### Tier 2: Strict (post-REINVENT → downstream analysis)

Applied after REINVENT generation for compounds proceeding to MD simulation.
Full 34-property assessment with standard literature thresholds.

**Usage:**
```bash
# Post-docking: keep more seeds
admet-filter --strictness relaxed --output seeds.json

# Post-REINVENT: rigorous filtering
admet-filter --strictness strict --output candidates.json
```
