# CMap/L1000 Drug Repurposing — Comprehensive Guide

> [!NOTE]
> This guide covers the complete theory, workflow, and biological application of the Connectivity Map (CMap) and L1000 platform for computational drug repurposing from SHAP-derived disease signatures. Based on: [clue.io/connectopedia](https://clue.io/connectopedia) | Subramanian et al., *Cell* 2017

---

## 1. What Is the Connectivity Map?

The **Connectivity Map (CMap)** is a resource developed at the **Broad Institute of MIT and Harvard** that enables discovery of functional connections between **drugs**, **genes**, and **diseases** through gene expression profiling.

**Core idea:** If a drug treatment produces transcriptional changes that are the **opposite** of a disease signature, that drug may reverse the disease state and serve as a therapeutic candidate.

### Key Capabilities

| Use Case | Description |
|---|---|
| **Drug Repurposing** | Find existing drugs that reverse a disease gene signature |
| **Mechanism of Action (MoA)** | Identify MoA of an uncharacterized compound by finding similar known perturbagens |
| **Target Discovery** | Connect genetic perturbations to drug effects to uncover novel targets |
| **Toxicity Prediction** | Identify compounds with signatures similar to known toxicants |

### Scale
- **>1.5 million** gene expression profiles
- **>25,000** perturbagens (small molecules, shRNA, cDNA, CRISPR)
- **~77** cell lines
- Generated using the **L1000** assay platform

### Key References

1. **Lamb J, et al.** (2006). The Connectivity Map: using gene-expression signatures to connect small molecules, genes, and disease. *Science*, 313(5795):1929–1935
2. **Subramanian A, et al.** (2017). A Next Generation Connectivity Map: L1000 Platform and the First 1,000,000 Profiles. *Cell*, 171(6):1437–1452
3. **Corsello SM, et al.** (2017). The Drug Repurposing Hub: a next-generation drug library. *Nature Medicine*, 23:405–408
4. **Corsello SM, et al.** (2020). Discovering the anticancer potential of non-oncology drugs. *Nature Cancer*, 1:235–248
5. **Pushpakom S, et al.** (2019). Drug repurposing: progress, challenges and recommendations. *Nature Reviews Drug Discovery*, 18:41–58
6. **He B, et al.** (2023). ASGARD: A Single-cell Guided pipeline for Drug repurposing. *Briefings in Bioinformatics*, 24(1):bbac593
7. **Kuleshov MV, et al.** (2016). Enrichr: comprehensive gene set enrichment analysis. *Nucleic Acids Research*, 44(W1):W90–W97

---

## 2. The L1000 Assay

The **L1000** is a high-throughput gene expression profiling technology based on Luminex bead technology. It directly measures **978 "landmark" genes** and computationally infers the expression of **~11,350 additional genes**.

### Gene Sets

| Gene Set | Count | Description |
|---|---|---|
| **Landmark (LM)** | 978 | Directly measured on L1000 beads |
| **Best Inferred (BING)** | ~10,174 | High-confidence inferred genes + landmarks |
| **Full Inferred** | ~12,328 | All inferred genes + landmarks |

### Data Processing Levels

| Level | Name | Description |
|---|---|---|
| **Level 1** | LXB | Raw fluorescence intensity |
| **Level 2** | GEX | Gene expression for 978 landmarks |
| **Level 3** | Q2NORM/INF | Quantile-normalized + inferred expression |
| **Level 4** | ZSPC/ZSVC | Differential expression signatures (robust z-scores vs. controls) |
| **Level 5** | MODZ | **Replicate-consensus signatures** — use this for queries |

> [!TIP]
> For CMap queries, always use **Level 5 (MODZ)** signatures — these are de-noised consensus profiles aggregated across biological replicates.

### Perturbation Types

| Code | Type | Description |
|---|---|---|
| **trt_cp** | Compound | Small molecule drug treatment |
| **trt_sh** | shRNA | Gene silencing via short hairpin RNA |
| **trt_oe** | Overexpression | Gene overexpression via cDNA/ORF |
| **trt_xpr** | CRISPR | Gene knockout via CRISPR |
| **ctl_vehicle** | Vehicle control | DMSO-treated control wells |

---

## 3. Understanding Connectivity Scores

### 3.1 Weighted Connectivity Score (WCS)

Based on a **weighted, signed two-sample Kolmogorov-Smirnov (KS) statistic**. Measures how enriched your query genes are at the top or bottom of each reference profile's ranked gene list.

### 3.2 Tau (τ) Score — The Key Metric

**Subramanian et al., Cell 2017**

Tau represents the **percentage of reference signatures** that have a less extreme connectivity score than the observed connection.

| Tau Score | Percentile | Interpretation | Action |
|---|---|---|---|
| **+90 to +100** | Top 1% similar | **Mimics** disease signature | Avoid (may worsen disease) |
| **+75 to +90** | Top 5% similar | Moderate mimic | Low priority |
| **-50 to +50** | Middle | No meaningful connection | Discard |
| **-75 to -90** | Top 5% reversed | Moderate reverser | Investigate |
| **-90 to -100** | **Top 1% reversed** | **Strongly reverses disease** | **Priority candidate** |

> [!IMPORTANT]
> **For drug repurposing:** Focus on perturbagens with **Tau ≤ -90**. These are compounds whose transcriptional effects are the **opposite** of your disease signature — meaning they may reverse the disease state.

### 3.3 Connectivity Across Cell Lines

- Scores computed **per cell line**
- **Summary Tau** (median across cell lines) assesses overall connection
- Cell-line-specific effects can reveal tissue-relevant biology (e.g., A549 for lung disease)

---

## 4. Multi-Level Drug Repurposing Pipeline

Modern drug repurposing integrates multiple computational layers for robust candidate identification:

```
Disease Signature (SHAP UP/DOWN) ─┬─► Level 1: Signature Reversal (CMap/L1000)
                                  ├─► Level 2: Pathway Enrichment (GSEA/GO/KEGG)
                                  ├─► Level 3: Network Pharmacology (PPI + Drug-Target)
                                  ├─► Level 4: Molecular Docking (optional)
                                  └─► Consensus Ranking → Prioritized Candidates
```

---

## 5. Level 1: Transcriptional Signature Reversal

### 5.1 Programmatic via GSEApy/Enrichr

**Kuleshov et al., Nucleic Acids Research 2016**

```python
import gseapy as gp

# Query L1000 perturbation databases
# Strategy: Find drugs that DOWNREGULATE disease-upregulated genes
results_up = gp.enrichr(
    gene_list=up_genes,
    gene_sets=[
        'L1000_Ligand_Pert_up',
        'L1000_Ligand_Pert_down',
        'Drug_Perturbations_from_GEO_up',
        'Drug_Perturbations_from_GEO_down',
        'LINCS_L1000_Chem_Pert_Consensus_Sigs',
    ],
    organism='Human',
    outdir='results/enrichr_up',
)

# Also query with DOWN genes
results_down = gp.enrichr(
    gene_list=down_genes,
    gene_sets=[
        'L1000_Ligand_Pert_up',
        'L1000_Ligand_Pert_down',
    ],
    organism='Human',
    outdir='results/enrichr_down',
)
```

### Available Gene Set Libraries

| Library | Description | Best For |
|---|---|---|
| `L1000_Ligand_Pert_up` | Genes upregulated by L1000 compounds | Primary query |
| `L1000_Ligand_Pert_down` | Genes downregulated by L1000 compounds | Primary query |
| `Drug_Perturbations_from_GEO_up` | Drug-induced upregulation (GEO) | Supplementary |
| `Drug_Perturbations_from_GEO_down` | Drug-induced downregulation (GEO) | Supplementary |
| `LINCS_L1000_Chem_Pert_Consensus_Sigs` | Consensus chemical perturbation sigs | High-confidence |
| `LINCS_L1000_Ligand_Perturbations_up` | Ligand perturbation signatures | Receptor-focused |

### 5.2 LINCS L1000 API (Direct)

```python
import requests

SIGCOM_URL = "https://maayanlab.cloud/sigcom-lincs"

# Submit signature query
payload = {
    "up_genes": up_genes[:150],    # Max 150 genes per list
    "down_genes": down_genes[:150],
    "direction": "reverse",
}
response = requests.post(f"{SIGCOM_URL}/api/v1/query", json=payload)
results = response.json()
```

### 5.3 Clue.io Manual Query — Step-by-Step

1. **Create account** at [clue.io](https://clue.io) (free for academic use)
2. Open **Query** app (`Ctrl+K` → Query)
3. Enter descriptive query name (e.g., "TB_granuloma_disease_signature")
4. Paste UP genes in the **UP** box (≤150 genes)
5. Paste DOWN genes in the **DOWN** box (≤150 genes)
6. System validates genes:
   - ✅ Valid, used: BING genes included
   - ⚠️ Valid, not used: Not in BING set
   - ❌ Invalid: Unrecognized
7. Click **Submit** → Wait ~5 minutes
8. Check **History** for results
9. Download results → Filter by **Tau ≤ -90**
10. Focus on **Level 5 (MODZ)** consensus signatures

### Signature Preparation Tips

- Use **HUGO gene symbols** (e.g., TP53, BRCA1) or Entrez Gene IDs
- Optimal: **50–150 genes per list** (top differentially expressed)
- Each list can have different lengths
- Can run with only UP or only DOWN list
- Genes must be in the BING gene set (~10,174 genes) to be used

### 5.4 Offline Fallback

```python
try:
    results = gp.enrichr(gene_list=genes, gene_sets=gene_set, organism='Human')
except Exception as e:
    print(f"API unreachable: {e}. Falling back to offline approach.")
    # Option 1: Load pre-downloaded gene set library
    # Available at: https://maayanlab.cloud/Enrichr/#libraries
    # Option 2: Use local CMap data from GEO (GSE92742)
    # Option 3: Submit to clue.io manually
```

---

## 6. Level 2: Pathway Enrichment Validation

**Subramanian et al., PNAS 2005**; **Kanehisa & Goto, NAR 2000**

Validate the disease signature and drug candidates through pathway analysis:

### Over-Representation Analysis (ORA)

```python
import gseapy as gp

# Enrichment of UP genes
up_enrichment = gp.enrichr(
    gene_list=up_genes,
    gene_sets=[
        'GO_Biological_Process_2023',
        'GO_Molecular_Function_2023',
        'KEGG_2021_Human',
        'Reactome_2022',
        'WikiPathways_2023_Human',
        'MSigDB_Hallmark_2020',
    ],
    organism='Human',
    outdir='results/enrichr_pathways_up',
)

# Enrichment of DOWN genes
down_enrichment = gp.enrichr(
    gene_list=down_genes,
    gene_sets=['GO_Biological_Process_2023', 'KEGG_2021_Human'],
    organism='Human',
    outdir='results/enrichr_pathways_down',
)
```

### Pre-Ranked GSEA (Full Gene Ranking)

Uses the complete SHAP-ranked gene list for more statistical power:

```python
# Use signed SHAP values as ranking metric
rnk = signature[['gene', 'mean_shap']].copy()
rnk.columns = ['Gene', 'Score']

gsea_results = gp.prerank(
    rnk=rnk,
    gene_sets=['MSigDB_Hallmark_2020', 'KEGG_2021_Human', 'Reactome_2022'],
    outdir='results/gsea_prerank',
    min_size=15,
    max_size=500,
    permutation_num=1000,
    seed=42,
)
```

---

## 7. Level 3: Network Pharmacology

**Barabási et al., Nature Reviews Genetics 2011**

### 7.1 Protein-Protein Interaction (PPI) Network

Build disease gene network from STRING database:

```python
import networkx as nx
import requests

def get_string_network(gene_list, species=9606, score_threshold=700):
    """Query STRING database for PPI network.
    Szklarczyk et al., Nucleic Acids Research 2023
    """
    url = "https://string-db.org/api/json/network"
    params = {
        "identifiers": "%0d".join(gene_list[:50]),  # API limit ~50 genes
        "species": species,
        "required_score": score_threshold,
        "network_type": "functional",
    }
    response = requests.get(url, params=params)
    if response.status_code != 200:
        print(f"STRING API error: {response.status_code}")
        return nx.Graph()
    
    interactions = response.json()
    G = nx.Graph()
    for edge in interactions:
        G.add_edge(edge["preferredName_A"], edge["preferredName_B"],
                    score=edge["score"])
    return G

# Build PPI network from top signature genes
ppi_network = get_string_network(top_genes["gene"].tolist()[:50])

# Network statistics
print(f"Nodes: {ppi_network.number_of_nodes()}")
print(f"Edges: {ppi_network.number_of_edges()}")
print(f"Density: {nx.density(ppi_network):.4f}")

# Hub genes (high centrality = key regulators)
degree_centrality = nx.degree_centrality(ppi_network)
betweenness_centrality = nx.betweenness_centrality(ppi_network)

hub_genes = sorted(degree_centrality.items(), key=lambda x: x[1], reverse=True)[:10]
print("\nTop 10 hub genes (by degree centrality):")
for gene, centrality in hub_genes:
    print(f"  {gene}: degree={centrality:.3f}, betweenness={betweenness_centrality.get(gene, 0):.3f}")
```

### 7.2 Drug-Target Interaction Databases

| Database | Content | URL | Reference |
|---|---|---|---|
| DrugBank | FDA-approved drugs + targets | drugbank.ca | Wishart et al., NAR 2018 |
| ChEMBL | Bioactivity data | ebi.ac.uk/chembl | Zdrazil et al., NAR 2024 |
| DGIdb | Drug-gene interactions | dgidb.org | Freshour et al., NAR 2021 |
| TTD | Therapeutic target database | db.idrblab.net/ttd | Zhou et al., NAR 2024 |
| STITCH | Chemical-protein interactions | stitch.embl.de | Szklarczyk et al., NAR 2016 |

```python
# Query DGIdb for known drug-gene interactions
def query_dgidb(gene_list):
    """Query DGIdb for known drug-gene interactions.
    Freshour et al., Nucleic Acids Research 2021
    """
    url = "https://dgidb.org/api/v2/interactions.json"
    params = {"genes": ",".join(gene_list[:50])}
    response = requests.get(url, params=params)
    if response.status_code != 200:
        return []
    
    data = response.json()
    interactions = []
    for match in data.get("matchedTerms", []):
        gene = match["searchTerm"]
        for interaction in match.get("interactions", []):
            interactions.append({
                "gene": gene,
                "drug": interaction.get("drugName", ""),
                "interaction_type": interaction.get("interactionTypes", ""),
                "source": interaction.get("sources", []),
            })
    return interactions

dgi_results = query_dgidb(top_genes["gene"].head(50).tolist())
print(f"Found {len(dgi_results)} known drug-gene interactions")
```

---

## 8. Level 4: Molecular Docking (Optional — Advanced)

For top candidates, validate binding to key disease targets:

```python
# Using AutoDock Vina for molecular docking
# Requires:
# 1. Protein structure (from PDB or AlphaFold2: alphafold.ebi.ac.uk)
# 2. Drug 3D structure (from PubChem or ChEMBL)
# 3. AutoDock Vina installed

# Command-line example:
# vina --receptor protein.pdbqt --ligand drug.pdbqt \
#      --center_x 10.0 --center_y 10.0 --center_z 10.0 \
#      --size_x 20 --size_y 20 --size_z 20 \
#      --exhaustiveness 32 --out result.pdbqt

# Binding affinity < -7.0 kcal/mol is considered strong binding
```

---

## 9. Single-Cell-Aware Drug Repurposing

**He et al., Briefings in Bioinformatics 2023** (ASGARD)

Traditional bulk approaches average across cell types. ASGARD uses cell-type-specific expression to identify drugs targeting specific cellular subpopulations:

```r
# ASGARD pipeline (R package)
# library(ASGARD)
# drug_results <- ASGARD_DrugRepurposing(
#     sce = single_cell_data,
#     cell_type = "macrophage",
#     disease_signature = signature_genes,
# )
```

---

## 10. Multi-Criteria Drug Candidate Prioritization

Score candidates across multiple evidence layers:

```python
def prioritize_candidates(enrichr_df, ppi_hub_set, dgi_drug_set):
    """
    Multi-criteria drug candidate scoring.
    
    Scoring rubric:
    - Signature reversal score (CMap/Enrichr):   0-40 points
    - Statistical significance (adj P-value):     0-20 points  
    - Target in PPI hub genes:                    0-15 points
    - Known drug-gene interaction (DGIdb):        0-15 points
    - FDA approval status:                        0-10 points
    """
    candidates = []
    max_score = enrichr_df['Combined Score'].max() if len(enrichr_df) > 0 else 1
    
    for _, row in enrichr_df.iterrows():
        score = 0
        drug_name = clean_drug_name(row['Term'])
        
        # Signature reversal strength (0-40)
        score += min(40, row['Combined Score'] / max_score * 40)
        
        # Statistical significance (0-20)
        adj_p = row.get('Adjusted P-value', 1.0)
        if adj_p < 0.001:
            score += 20
        elif adj_p < 0.01:
            score += 15
        elif adj_p < 0.05:
            score += 10
        
        # Target overlaps with PPI hub genes (0-15)
        overlap_genes = set(row.get('Genes', '').split(';')) & ppi_hub_set
        if len(overlap_genes) > 0:
            score += min(15, len(overlap_genes) * 5)
        
        # Known interaction in DGIdb (0-15)
        if drug_name.lower() in dgi_drug_set:
            score += 15
        
        candidates.append({
            'drug': drug_name,
            'raw_term': row['Term'],
            'combined_score': row['Combined Score'],
            'adj_pvalue': adj_p,
            'overlap_genes': ';'.join(overlap_genes) if overlap_genes else '',
            'priority_score': score,
        })
    
    candidates = sorted(candidates, key=lambda x: x['priority_score'], reverse=True)
    return pd.DataFrame(candidates)
```

---

## 11. Drug Name Cleaning

L1000/Enrichr terms include cell line and direction suffixes:

```python
import re

def clean_drug_name(name):
    """Extract compound name from L1000/Enrichr term."""
    # Remove cell line suffixes
    name = re.sub(r'\s+(HL60|PC3|MCF7|A375|A549|HCC515|HEPG2|HT29|VCAP)\s+(UP|DOWN)$',
                  '', name, flags=re.IGNORECASE)
    # Remove database IDs
    name = re.sub(r'\s+(DB\d+|TTD\s*\d+|CTD\s*\d+)\s*', ' ', name)
    # Remove GEO sample info
    name = re.sub(r'\s+(human|rat|mouse)\s+GSE\d+\s+sample\s+\d+', '', name)
    # Remove trailing numbers
    name = re.sub(r'[\s-]+\d{3,}$', '', name)
    return ' '.join(name.split()).strip().capitalize()
```

---

## 12. Touchstone — The Reference Dataset

### Composition
- **~2,400 compounds** + **~3,700 genetic perturbagens**
- Core set of **9 cell lines**: A375, A549, HA1E, HCC515, HEPG2, HT29, MCF7, PC3, VCAP
- Well-annotated with MoA, target, and clinical information

### Cell Line Relevance for Disease Studies

| Cell Line | Origin | Relevant For |
|---|---|---|
| **A549** | Lung adenocarcinoma | TB (lung), COPD, lung cancer |
| **MCF7** | Breast cancer | Most commonly profiled (largest coverage) |
| **HCC515** | Non-small cell lung cancer | Respiratory diseases |
| **HEPG2** | Hepatocellular carcinoma | Liver diseases, drug metabolism |
| **HT29** | Colorectal adenocarcinoma | GI diseases |
| **PC3** | Prostate cancer | Prostate cancer |

### Query Modes

| Mode | When to Use |
|---|---|
| **Unmatched (recommended)** | Compares against all Touchstone signatures regardless of cell type |
| **Matched** | Restricts to matching cell types (only if query from core 9 cell lines) |

---

## 13. Programmatic Access

### cmapPy (Python)

```python
# pip install cmapPy
from cmapPy.pandasGEXpress import parse

gctoo = parse.parse("my_data.gctx")
print(gctoo.data_df.shape)          # (genes, samples)
print(gctoo.row_metadata_df.head()) # gene annotations
print(gctoo.col_metadata_df.head()) # sample annotations
```

### Data Access

CMap/LINCS L1000 data available from:
1. **[clue.io/data](https://clue.io/data)** — Primary portal
2. **GEO** — GSE92742 (Phase I), GSE70138 (Phase II)
3. **LINCS Data Portal** — [lincsproject.org](http://lincsproject.org)

### File Formats

| Format | Description |
|---|---|
| **GCT** | Tab-delimited annotated matrix |
| **GCTX** | HDF5-based binary format (large datasets) |
| **GMT** | Gene Matrix Transposed — gene set format |

---

## 14. Validation & Experimental Follow-Up

### Computational Validation

1. **Cross-dataset validation:** Test signature on independent cohorts (different GEO accessions)
2. **Leave-one-out analysis:** Remove each top gene and re-run pipeline to test stability
3. **Random signature control:** Compare drug repurposing results to 1,000 random gene sets (permutation test)
4. **Cross-reference CMap with literature:** PubMed search for "[drug_name] AND [disease]"

### Experimental Validation Hierarchy

1. **In vitro assays:** Treat relevant cell lines with top candidates
2. **Dose-response curves:** IC50/EC50 determination
3. **Transcriptomic validation:** RNA-seq after drug treatment to confirm signature reversal
4. **Animal models:** Pre-clinical efficacy and safety
5. **Clinical trials:** Translational validation

---

## 15. CMap Glossary

| Term | Definition |
|---|---|
| **Perturbagen** | Any agent used to perturb cells (compounds, shRNA, cDNA, CRISPR) |
| **Signature** | Vector of differential expression z-scores induced by a perturbation |
| **Landmark genes** | The 978 genes directly measured by L1000 |
| **BING** | Best INferred Genes (~10,174 well-inferred + 978 landmarks) |
| **Connectivity** | Quantitative similarity between two expression signatures |
| **Tau (τ)** | Normalized connectivity score (-100 to +100) |
| **TAS** | Transcriptional Activity Score (perturbation strength + replicate concordance) |
| **MODZ** | Moderated Z-score — replicate aggregation method |
| **Introspect** | Querying a dataset against itself |
| **Robust z-score** | Z-score using median and MAD (outlier-resistant) |

---

## 16. Key References

1. Lamb J, et al. (2006). The Connectivity Map. *Science*, 313(5795):1929–1935
2. Subramanian A, et al. (2017). A Next Generation Connectivity Map. *Cell*, 171(6):1437–1452
3. Corsello SM, et al. (2017). The Drug Repurposing Hub. *Nat Med*, 23:405–408
4. Corsello SM, et al. (2020). Anticancer potential of non-oncology drugs. *Nat Cancer*, 1:235–248
5. Pushpakom S, et al. (2019). Drug repurposing. *Nat Rev Drug Discov*, 18:41–58
6. He B, et al. (2023). ASGARD for single-cell drug repurposing. *Brief Bioinform*, 24(1):bbac593
7. Kuleshov MV, et al. (2016). Enrichr. *NAR*, 44(W1):W90–W97
8. Barabási A-L, et al. (2011). Network medicine. *Nat Rev Genet*, 12:56–68
9. Szklarczyk D, et al. (2023). STRING database 2023. *NAR*, 51(D1):D483–D489
10. Wishart DS, et al. (2018). DrugBank 5.0. *NAR*, 46(D1):D1074–D1082
11. Freshour SL, et al. (2021). Integration of the DGIdb. *NAR*, 49(D1):D1144–D1151
12. Subramanian A, et al. (2005). Gene set enrichment analysis. *PNAS*, 102(43):15545–15550

## Support & Resources

| Resource | Link |
|---|---|
| **Clue Platform** | [clue.io](https://clue.io) |
| **Connectopedia** | [clue.io/connectopedia](https://clue.io/connectopedia) |
| **cmapPy (Python)** | [github.com/cmap/cmapPy](https://github.com/cmap/cmapPy) |
| **cmapR (R)** | [github.com/cmap/cmapR](https://github.com/cmap/cmapR) |
| **Morpheus** | [software.broadinstitute.org/morpheus](https://software.broadinstitute.org/morpheus) |
| **Email Support** | clue@broadinstitute.org |
