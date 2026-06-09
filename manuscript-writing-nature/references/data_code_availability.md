# Data & Code Availability — FAIR Principles

## Overview

Nature mandates Data Availability and Code Availability statements for all manuscripts. These must follow the **FAIR principles**: Findable, Accessible, Interoperable, and Reusable.

---

## Data Availability Statement

### Required Elements
1. **Where** the data can be accessed (repository name + URL)
2. **Identifier** (accession number, DOI, or persistent URL)
3. **Access conditions** (open, restricted, or upon request — with justification)
4. **Source data** for all figures (provided with the paper or in repository)

### Approved Repositories

| Data Type | Repository | URL |
|---|---|---|
| Sequencing (RNA-seq, scRNA-seq) | GEO / SRA | ncbi.nlm.nih.gov/geo |
| Sequencing (European) | ArrayExpress / ENA | ebi.ac.uk/arrayexpress |
| Proteomics | PRIDE | ebi.ac.uk/pride |
| Metabolomics | MetaboLights | ebi.ac.uk/metabolights |
| Structural biology | PDB / EMDB | rcsb.org |
| General purpose | Zenodo | zenodo.org |
| General purpose | Figshare | figshare.com |
| General purpose | Dryad | datadryad.org |
| Genomic variants | ClinVar / dbSNP | ncbi.nlm.nih.gov/clinvar |
| Clinical trials | ClinicalTrials.gov | clinicaltrials.gov |

### Template Statements

**Open data (preferred)**:
> The single-cell RNA-sequencing data generated in this study have been deposited in the Gene Expression Omnibus (GEO) under accession number GSE178837. Processed count matrices and cell-level metadata are available at Zenodo (https://doi.org/10.5281/zenodo.XXXXXXX). Source data for Figs. 1–6 and Extended Data Figs. 1–10 are provided with this paper.

**Previously published data**:
> This study used publicly available single-cell RNA-seq data from GSE178837 (ref. X). No new data were generated. All analysis outputs are available at [repository URL].

**Restricted data (must justify)**:
> The clinical data used in this study contain protected health information and cannot be shared publicly. De-identified summary statistics are available at [repository]. Individual-level data are available from [contact] upon reasonable request and execution of a data use agreement, in compliance with [IRB/ethics board].

**NOT acceptable** (Nature will reject):
> ❌ "Data are available upon request."
> ❌ "Data are available from the corresponding author."
> ❌ No statement at all.

---

## Code Availability Statement

### Required Elements
1. **Repository URL** (GitHub, GitLab, Bitbucket)
2. **Archived version with DOI** (Zenodo, Code Ocean, Software Heritage)
3. **Version/commit** used in the manuscript
4. **License** (MIT, Apache-2.0, GPL-3.0, etc.)
5. **Dependencies** and installation instructions

### Best Practices
- Archive the code on **Zenodo** (not just GitHub — GitHub repos can be deleted)
- Create a **release tag** matching the manuscript version
- Include a **README** with installation and usage instructions
- Include a **requirements.txt** or **environment.yml**
- Include **test data** or a small example that reproduces key results
- Use a **DOI** for the archived version

### How to Archive GitHub on Zenodo
1. Go to https://zenodo.org → Log in with GitHub
2. Enable the repository in Zenodo settings
3. Create a GitHub release (tag: `v1.0.0-manuscript`)
4. Zenodo automatically creates an archived copy with DOI
5. Use the Zenodo DOI in your manuscript

### Template Statements

**Custom analysis code**:
> Custom code for the ensemble machine learning pipeline, SHAP analysis, and drug repurposing workflow is available at https://github.com/username/project (ref. X) and archived at Zenodo (https://doi.org/10.5281/zenodo.XXXXXXX). The analysis was performed using Python 3.11 with scikit-learn v1.3, XGBoost v1.7, LightGBM v4.0, SHAP v0.42, and Scanpy v1.9.

**Standard tools only**:
> No custom code was developed for this study. All analyses used publicly available software packages as described in the Methods section.

---

## FAIR Principles Checklist

### Findable
```
□ Data/code has a globally unique persistent identifier (DOI, accession number)
□ Data/code is described with rich metadata
□ Metadata clearly includes the identifier of the data it describes
□ Data/code is registered or indexed in a searchable resource
```

### Accessible
```
□ Data/code is retrievable by its identifier using a standardized protocol
□ The protocol is open, free, and universally implementable
□ Metadata is accessible even when the data is not (embargo, restricted)
□ Authentication/authorization procedures are clearly documented if needed
```

### Interoperable
```
□ Data uses a formal, accessible, shared, and broadly applicable language
□ Data uses vocabularies that follow FAIR principles (ontologies, controlled vocabularies)
□ Data includes qualified references to other (meta)data
□ File formats are standard and open (CSV, FASTQ, BAM, HDF5, not proprietary)
```

### Reusable
```
□ Data is released with a clear and accessible data usage license
□ Data is associated with detailed provenance (how it was generated)
□ Data meets domain-relevant community standards (MIAME, MINSEQE, etc.)
□ Code includes documentation, tests, and examples
```

---

## Checklist Before Submission

```
□ Data deposited in appropriate repository
□ Accession numbers / DOIs obtained
□ Data Availability statement written with specific identifiers
□ Code archived on Zenodo (or equivalent) with DOI
□ Code repository includes README, LICENSE, requirements.txt
□ Code Availability statement written with URLs and versions
□ Source data for all figures provided
□ All software versions listed in Methods
□ FAIR principles met for all shared resources
```
