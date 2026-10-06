# Example REINVENT 4 Configurations

These files were produced by `ag-compound-synthesis generate-config` and
each one parses as valid TOML. They are starting points: adjust the paths,
the device and the step counts for your run, or regenerate them with the
command shown at the top of each file's header comment.

| File | Mode | Generator | Purpose |
|---|---|---|---|
| `denovo_rl.toml` | staged_learning | Reinvent | Unconstrained de novo design against drug-like properties |
| `scaffold_hopping.toml` | staged_learning | Mol2Mol | Hop to new scaffolds from a known lead |
| `transfer_learning.toml` | transfer_learning | Reinvent | Focus the prior onto a set of known actives |
| `anti_tb_rl.toml` | staged_learning | Reinvent | Anti-tubercular property space, seeded from screening hits |
| `scoring_only.toml` | scoring | - | Score an existing set without generating |

Prefer regenerating over editing by hand: the generator validates the
combination, resolves the prior path, refuses a generator whose seeds are
missing, and writes a `.meta.json` recording the resolved settings.
