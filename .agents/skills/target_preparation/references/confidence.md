# Structure Confidence Reference

Deep detail for `assess-structure`. The SKILL.md carries what you need to
run it; read this when a verdict is surprising or needs defending.

## Three B-factor scales, told apart

The B-factor column means different things in different files, and guessing
wrong produces a confident, wrong verdict in either direction. The `scale`
field in the output says which was detected.

| `scale` | Source | Convention |
|---|---|---|
| `plddt_0_100` | AlphaFold DB, ColabFold | Higher is better, bounded at 100 |
| `plddt_0_1` | The **ESMFold API** | The same quantity as a fraction; rescaled by 100 on read |
| `bfactor` | Experimental structure | Atomic displacement parameters: unbounded, **lower is better** |

### Why the fractional case matters

The ESMFold API writes pLDDT between 0 and 1. Measured against the live
endpoint on a 122-residue fold: B-factors from 0.40 to 0.94, mean 0.842.
That is an 84 pLDDT model — a good one. Read as a 0–100 score it becomes
0.84, which any sane banding calls "very low confidence, not recommended for
docking". Detection is therefore not a nicety: without it a usable model is
rejected outright.

Detection rule: every value in [0, 1] across at least 5 residues means a
fraction. A maximum above 100 means crystallographic. Otherwise a mean at or
above 20 means a 0–100 pLDDT scale, and below that, crystallographic — a
well-ordered crystal structure clusters low, a usable predicted model does
not.

## pLDDT bands

Following the AlphaFold DB convention (Jumper et al. 2021; Varadi et al.
2022):

| Band | Range | Interpretation |
|---|---|---|
| Very high | >= 90 | Backbone and side chains both reliable |
| Confident | 70–90 | Backbone reliable; side chains less so |
| Low | 50–70 | Treat with caution; often a flexible region |
| Very low | < 50 | Should not be interpreted as structure |

## How the verdict is reached

```
scale is bfactor                      -> NOT_APPLICABLE
mean < 50  OR  >50% of residues < 50  -> STOP
mean < 70  OR  a disordered run of 5+ residues
           OR  >30% of residues in 50-70                -> REFINE
otherwise                                               -> PROCEED
```

A `STOP` verdict lists three routes: re-predict with an MSA-based method,
search the PDB for a homologue, or fall back to ligand-based screening.

A `REFINE` verdict names the disordered regions to truncate and suggests
energy minimisation before docking.

## Low-confidence regions

Reported as contiguous runs below pLDDT 50, per chain. A run is broken by a
chain change **or** a gap in residue numbering, so a region is always
genuinely contiguous in sequence. Crystal structures have unmodelled
stretches — residue 50 can be followed by residue 80 — and treating those as
adjacent would report a region that does not exist.

Keep the docking box away from these regions. The `actions` field says so
explicitly when any are found.

## Multi-chain structures

Confidence is keyed on `(chain, residue number, insertion code)`. Keying on
the residue number alone — as an earlier version of this skill did — means
residue 42 of chain B overwrites residue 42 of chain A, so a dimer reports
the confidence of roughly half its residues. Per-chain means are in
`per_chain`.

## Experimental structures

When `scale` is `bfactor`, judge the structure by:

-   **Resolution** — below about 2.5 Å for reliable side-chain placement.
-   **R-free** — and the gap between R and R-free.
-   **Occupancy and B-factors of the binding-site residues specifically**, not
    the global mean.
-   **Completeness** — `pdbqt.warnings` from `prepare-receptor` names
    residues with missing atoms, and flags any inside the docking box.

## Output size

`assess-structure` writes the full report. Other commands that embed it
(`prepare-receptor`) carry only the summary by default, because the full
form is roughly five times larger and is mostly per-region detail. Pass
`--detail full` to embed everything.
