# REINVENT4 Prior Models Reference

## Available Models

All priors are available on [Zenodo](https://doi.org/10.5281/zenodo.15641296).

| Key | Filename | Architecture | Use Case |
|---|---|---|---|
| `.reinvent` | `reinvent.prior` | GRU-RNN | De novo generation (no seeds needed) |
| `.libinvent` | `libinvent.prior` | Transformer | R-group replacement on scaffolds |
| `.linkinvent` | `linkinvent.prior` | Transformer | Linker design between fragments |
| `.m2m_scaffold` | `mol2mol_scaffold.prior` | Transformer | Scaffold hopping |
| `.m2m_scaffold_generic` | `mol2mol_scaffold_generic.prior` | Transformer | Generic scaffold optimization |
| `.m2m_high` | `mol2mol_high_similarity.prior` | Transformer | High-similarity analogs |
| `.m2m_medium` | `mol2mol_medium_similarity.prior` | Transformer | Medium-similarity analogs |
| `.m2m_mmp` | `mol2mol_mmp.prior` | Transformer | Matched molecular pair transforms |
| `.pepinvent` | `pepinvent.prior` | Transformer | Peptide design |

## Choosing a Generator

```
Want to design from scratch?
  └─ Use Reinvent (.reinvent) — no seeds needed

Want to optimize a known drug?
  └─ How similar should results be?
       ├─ Very similar → Mol2Mol (.m2m_high)
       ├─ Moderately similar → Mol2Mol (.m2m_medium)
       └─ Different scaffold → Mol2Mol (.m2m_scaffold_generic)

Want to keep the core scaffold, change R-groups?
  └─ Use LibInvent (.libinvent)

Want to link two fragments?
  └─ Use LinkInvent (.linkinvent)

Want to design peptides?
  └─ Use Pepinvent (.pepinvent)
```

## Download

```bash
# Via compound_synthesis.py
python compound_synthesis.py download-priors \
    --reinvent-dir E:\ANTIGRAVITY_WORKSHOP\REINVENT4 \
    --models reinvent libinvent \
    --output download_report.json

# Manual (curl)
# Visit https://doi.org/10.5281/zenodo.15641296
# Download .prior files into REINVENT4/priors/
```

## Custom Prior Location

Set the `REINVENT_PRIOR_BASE` environment variable to use priors from a
different directory:

```bash
export REINVENT_PRIOR_BASE=/path/to/my/priors
```

Or use dot notation in TOML configs to reference built-in registry:
```toml
prior_file = ".reinvent"  # resolved at runtime
```

## Seed File Formats

| Generator | File Format | Example |
|---|---|---|
| Reinvent | Not required | — |
| LibInvent | 1 scaffold/line, `*` marks attachment | `c1ccc(*)cc1N(*)` |
| LinkInvent | 2 warheads/line, `\|` separator | `c1ccccc1N\|Nc1ccccc1` |
| Mol2Mol | 1 SMILES/line | `CC(=O)Oc1ccccc1C(=O)O` |
| Pepinvent | 1 peptide SMILES/line | `CC(N)C(=O)NC(CC(=O)O)C(=O)O` |
