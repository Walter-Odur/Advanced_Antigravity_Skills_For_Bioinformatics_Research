# Exemplary Nature Articles — Study These as Models

## Purpose

These are real, high-impact articles published in Nature that exemplify the writing standards, narrative structure, and presentation quality that every manuscript should aspire to. Study their abstracts, introductions, figure designs, and discussion structures.

---

## Computational Biology / Machine Learning in Biology

### 1. AlphaFold: Protein Structure Prediction
**Jumper, J. et al. Highly accurate protein structure prediction with AlphaFold. Nature 596, 583–589 (2021).**
- DOI: 10.1038/s41586-021-03819-2
- Citations: 25,000+
- **Why study this paper:**
  - Masterclass in making computational methods accessible to broad audience
  - Abstract opens with "Proteins are essential to life..." — hooks every scientist
  - Figures are self-explanatory: each panel tells a complete story
  - Methods are exhaustively detailed for reproducibility
  - Discussion is concise but places work in transformative context

### 2. The Connectivity Map (CMap L1000)
**Subramanian, A. et al. A next generation connectivity map: L1000 platform and the first 1,000,000 profiles. Cell 171, 1437–1452 (2017).**
- DOI: 10.1016/j.cell.2017.10.049
- Citations: 3,500+
- **Why study this paper:**
  - Resource paper that clearly explains both the platform AND the science
  - Excellent figures that visualize high-dimensional data intuitively
  - Methods section is a gold standard for computational reproducibility
  - Validation through multiple independent approaches

### 3. Single-Cell Best Practices
**Heumos, L. et al. Best practices for single-cell analysis across modalities. Nat. Rev. Genet. 24, 550–572 (2023).**
- DOI: 10.1038/s41576-023-00586-w
- Citations: 500+
- **Why study this paper:**
  - Comprehensive review structure — systematic and thorough
  - Figures are conceptual diagrams that explain workflows
  - Accessible to newcomers while being useful to experts
  - Clear recommendations backed by benchmarking evidence

### 4. SHAP: Interpretable Machine Learning
**Lundberg, S. M. et al. From local explanations to global understanding with explainable AI for trees. Nat. Mach. Intell. 2, 56–67 (2020).**
- DOI: 10.1038/s42256-019-0138-9
- Citations: 4,000+
- **Why study this paper:**
  - Makes complex mathematics accessible to biologists
  - Figures demonstrate the method on real biomedical data
  - Bridges CS methodology with clinical application
  - Clear problem → solution → validation structure

### 5. Scanpy: Single-Cell Analysis in Python
**Wolf, F. A., Angerer, P. & Theis, F. J. SCANPY: large-scale single-cell gene expression data analysis. Genome Biol. 19, 15 (2018).**
- DOI: 10.1186/s13059-017-1382-0
- Citations: 6,000+
- **Why study this paper:**
  - Software paper that balances technical detail with accessibility
  - Benchmarks against existing tools with fair comparisons
  - Figures show real biological results, not just benchmarks
  - Methods are reproducible with provided code

---

## Drug Repurposing / Therapeutics

### 6. Host-Directed Therapy for TB
**Kaufmann, S. H. E. et al. Host-directed therapies for bacterial and viral infections. Nat. Rev. Drug Discov. 17, 35–56 (2018).**
- DOI: 10.1038/nrd.2017.162
- Citations: 1,000+
- **Why study this paper:**
  - Comprehensive review of host-directed therapeutic strategies
  - Excellent integration of immunology and pharmacology
  - Tables summarize drug candidates with clinical status
  - Discussion of challenges and future directions is balanced

### 7. Computational Drug Repurposing
**Pushpakom, S. et al. Drug repurposing: progress, challenges and recommendations. Nat. Rev. Drug Discov. 18, 41–58 (2019).**
- DOI: 10.1038/nrd.2018.168
- Citations: 3,500+
- **Why study this paper:**
  - Defines the field of drug repurposing clearly
  - Covers computational AND experimental approaches
  - Honest discussion of failure modes and limitations
  - Practical recommendations for the field

---

## Immunology / Tuberculosis

### 8. TB Granuloma Biology
**Gideon, H. P. et al. Multimodal profiling of lung granulomas in macaques reveals cellular correlates of tuberculosis control. Immunity 55, 827–846 (2022).**
- DOI: 10.1016/j.immuni.2022.04.004
- Citations: 200+
- **Why study this paper:**
  - Multi-modal data integration (scRNA-seq + spatial + flow cytometry)
  - Figures combine multiple data modalities in single panels
  - Biological interpretations are grounded in immunological literature
  - Statistical rigor with appropriate corrections

### 9. Single-Cell Atlas of the Immune System
**Dominguez Conde, C. et al. Cross-tissue immune cell analysis reveals tissue-specific features in humans. Science 376, eabl5197 (2022).**
- DOI: 10.1126/science.abl5197
- Citations: 1,000+
- **Why study this paper:**
  - Large-scale single-cell atlas with clear biological narrative
  - Figures are publication-quality with consistent color schemes
  - Methods are reproducible with public data and code
  - Integration across tissues with meaningful biological conclusions

---

## What to Learn from Each Paper

### Abstract Writing
Study the AlphaFold abstract — it opens with universal significance ("Proteins are essential to life"), narrows to the specific problem (structure prediction), states the approach (deep learning), and delivers the punchline (atomic-level accuracy). All in ~150 words.

### Introduction Structure
Study the CMap paper — the introduction follows a perfect funnel: broad significance of perturbational biology → specific gap in profiling scale → their solution (L1000 platform) → preview of what the paper delivers.

### Figure Design
Study Gideon et al. — multi-panel figures that combine different data types (UMAP, heatmaps, bar charts) with consistent colors and clear legends. Each figure tells a complete story.

### Discussion Quality
Study Pushpakom et al. — the discussion honestly acknowledges challenges in drug repurposing (high attrition rates, intellectual property issues) while making a constructive case for the field's future.

### Methods Reproducibility
Study the Scanpy paper — every parameter, every software version, every dataset accession number is documented. A reader can reproduce the entire analysis.

---

## Style Patterns to Emulate

### Opening Sentences (Hooks)

From these exemplary papers:

| Paper | Opening Sentence |
|---|---|
| AlphaFold | "Proteins are essential to life, and understanding their structure can facilitate a mechanistic understanding of their function." |
| CMap | "A fundamental challenge in biomedicine is to understand how diseases arise and how they may be treated." |
| Heumos et al. | "Single-cell technologies have provided unprecedented insights into cell types and states." |
| Kaufmann et al. | "Infectious diseases remain a major threat to human health worldwide." |
| Pushpakom et al. | "Drug repurposing — the process of finding new uses for existing drugs — has gained considerable momentum." |

**Pattern**: Start with a universal truth that every scientist would agree with, then narrow.

### "Here we..." Statements

| Paper | Contribution Statement |
|---|---|
| AlphaFold | "Here we present AlphaFold, a system for predicting protein structure..." |
| CMap | "Here we describe the generation of the next generation Connectivity Map..." |
| Scanpy | "Here we present Scanpy, a scalable toolkit for analyzing single-cell..." |

**Pattern**: Always starts with "Here we [verb]" — present, describe, demonstrate, show, report.

### Results Paragraph Structure

From AlphaFold (typical pattern):
```
[Setup sentence connecting to previous paragraph].
We [action verb] and found that [main result] (Fig. X).
[Quantitative detail supporting the result].
[Comparison to prior work or expected performance].
[Brief interpretation or significance of the result].
```

### Discussion Limitation Acknowledgement

From Gideon et al. (honest and constructive):
```
"Several limitations of our study should be noted. First, [limitation 1
with context for why it exists]. Second, [limitation 2]. Despite these
limitations, [constructive framing of what the work achieves despite
the limitations]. Future studies incorporating [proposed solutions]
could address these constraints."
```
