# External API Endpoints Reference

## RCSB Protein Data Bank (PDB)

Base URL: `https://files.rcsb.org`

- **Download PDB**: `GET /download/{PDB_ID}.pdb`
- **Download CIF**: `GET /download/{PDB_ID}.cif`

Metadata API: `https://data.rcsb.org/rest/v1/core`

- **Entry**: `GET /entry/{PDB_ID}`
- **Polymer entity**: `GET /polymer_entity/{PDB_ID}/{ENTITY_ID}`

## AlphaFold Database

Base URL: `https://alphafold.ebi.ac.uk`

- **PDB model**: `GET /files/AF-{UNIPROT_ID}-F1-model_v4.pdb`
- **CIF model**: `GET /files/AF-{UNIPROT_ID}-F1-model_v4.cif`
- **PAE JSON**: `GET /files/AF-{UNIPROT_ID}-F1-predicted_aligned_error_v4.json`
- **API metadata**: `GET /api/prediction/{UNIPROT_ID}`

## ChEMBL REST API

Base URL: `https://www.ebi.ac.uk/chembl/api/data`

### Common Endpoints

| Endpoint | Path | Searchable |
|----------|------|------------|
| Activity | `/activity.json` | Yes |
| Target | `/target.json` | Yes |
| Molecule | `/molecule.json` | Yes |
| Assay | `/assay.json` | Yes |

### Filter Operators (Django-style)

| Operator | Description | Example |
|----------|-------------|---------|
| (none) | Exact match | `target_chembl_id=CHEMBL203` |
| `__gte` | Greater or equal | `pchembl_value__gte=7` |
| `__lte` | Less or equal | `mw_freebase__lte=500` |
| `__in` | Value in list | `assay_type__in=B,F` |
| `__isnull` | Null check | `pchembl_value__isnull=false` |
| `__range` | Range | `mw_freebase__range=200,500` |

### Search

`GET /target/search.json?q=EGFR&limit=5`

### Pagination

All list endpoints return paginated results with `limit` and `offset` params.
Response includes `page_meta` with `total_count`.

## ZINC Database

Base URL: `https://zinc15.docking.org`

- **Substances**: `GET /substances.json?count=N`
- **Subsets**: drug-like, lead-like, fragment-like
- Filters: `mwt__lte`, `logp__lte`, `reactive`

## Useful jq Patterns

```bash
# Extract SMILES from ChEMBL activity results
jq '.compounds[].canonical_smiles' results.json

# Get top 5 compounds by pChEMBL
jq '.compounds | sort_by(-.pchembl_value) | .[0:5]' results.json

# List ADMET-passing compounds
jq '.compounds[] | select(.overall_pass==true) | .name' admet.json

# Get docking scores
jq '.ranked_results[] | {name, score: .score_kcal_mol}' dock.json
```
