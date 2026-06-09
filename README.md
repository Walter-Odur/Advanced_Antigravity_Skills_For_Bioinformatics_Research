# 🧬 Advanced Antigravity Skills for Bioinformatics Research

[![Skills](https://img.shields.io/badge/Skills-2-green.svg)](#available-skills)
[![Citations](https://img.shields.io/badge/Citations-70%2B-orange.svg)](#references)

A curated collection of **advanced [Antigravity IDE](https://github.com/google-deepmind/antigravity) skills** purpose-built for computational biology, bioinformatics, and translational research. Each skill is a comprehensive, literature-backed pipeline that encodes best practices from high-impact peer-reviewed journals (*Nature*, *Cell*, *Nature Methods*, *Genome Biology*, *Nature Machine Intelligence*).

> **More skills coming soon.** This repository is actively maintained and will grow with new bioinformatics skills covering spatial transcriptomics, multi-omics integration, CRISPR screen analysis, and more.

---

## 🎯 What Are Antigravity Skills?

[Antigravity](https://github.com/google-deepmind/antigravity) is an agentic AI coding assistant. **Skills** are structured instruction packages that extend its capabilities for specialized domains. When you install a skill, Antigravity gains deep knowledge of:

- Domain-specific best practices and methodologies
- Complete code templates and reference implementations
- Literature-backed parameter choices with citations
- Common pitfalls and troubleshooting guides

Think of skills as **expert knowledge modules** that transform a general-purpose AI into a domain specialist.

---

## 📦 Available Skills

### 1. `scrna-drug-repurposing`

**End-to-end single-cell RNA-seq ensemble ML pipeline for disease signature extraction and computational drug repurposing.**

| Metric | Value |
|---|---|
| **Lines of code** | 5,398 |
| **Files** | 16 |
| **Peer-reviewed citations** | 60+ |
| **Journals** | *Nature*, *Cell*, *Nature Methods*, *Genome Biology*, *Nature Machine Intelligence*, *Nature Reviews Drug Discovery* |

#### Pipeline Overview

```
scRNA-seq Data → QC & Preprocessing → Ensemble ML Training → SHAP Signature → Drug Repurposing → Manuscript Figures
```

#### What It Does

| Step | Description | Key Tools |
|---|---|---|
| **1. QC & Preprocessing** | Tiered QC with ambient RNA removal, doublet detection, normalization, batch correction | Scanpy, CellBender, Scrublet, Harmony, scVI |
| **2. ML Data Preparation** | Binary label mapping, stratified train/test splitting, feature extraction | scikit-learn, NumPy |
| **3. Ensemble ML Training** | Random Forest, XGBoost, LightGBM, Stacking Classifier with stratified CV | scikit-learn, XGBoost, LightGBM |
| **4. Model Improvement** | Hyperparameter tuning, learning curves, baseline comparison | Optuna, scikit-learn |
| **5. SHAP Analysis** | Interventional SHAP (Pearl's do-calculus), disease signature extraction | SHAP, TreeExplainer |
| **6. Drug Repurposing** | CMap/L1000 signature reversal, pathway enrichment, network pharmacology | GSEApy, Enrichr, STRING |
| **7. Manuscript Figures** | Publication-quality multi-panel figures (Nature/Cell/Science formatting) | Matplotlib, Seaborn |

#### Reference Documentation

| Document | Lines | Coverage |
|---|---|---|
| `qc_preprocessing.md` | 574 | Complete Scanpy user guide: AnnData, data loading, tiered QC, normalization, batch correction, clustering, annotation, trajectory, DE analysis |
| `shap_interpretation.md` | 668 | Interventional SHAP, interaction values, bootstrap stability, SHAP clustering, pathway enrichment validation |
| `drug_repurposing.md` | 486 | 4-level evidence pipeline: CMap signature reversal, GSEA pathway enrichment, PPI network pharmacology, molecular docking |
| `ensemble_ml.md` | 392 | Model configs, class imbalance, nested CV, Optuna tuning, probability calibration, biology-specific ML patterns |
| `manuscript_figures.md` | 260 | Nature/Cell/Science column widths, figure recipes, statistical annotations, color-blind-safe palettes |
| `pipeline_architecture.md` | 130 | Compute estimates, GPU requirements, parallelization strategy, dependency graph |

#### Key Methodological Choices

- **Interventional SHAP** over observational — uses Pearl's do-calculus to handle correlated gene expression features (Lundberg et al., *Nat Mach Intell* 2020)
- **Tau ≤ -90** threshold for CMap drug candidates — top 1% signature reversers (Subramanian et al., *Cell* 2017)
- **MCC** over accuracy/F1 for imbalanced data — considers all four confusion matrix quadrants (Chicco & Jurman, *BMC Genomics* 2020)
- **Nested cross-validation** for unbiased performance estimation when tuning (Cawley & Talbot, *JMLR* 2010)
- **Tissue-specific MT% thresholds** — 5% for PBMCs, 20% for lung/liver (Heumos et al., *Nat Rev Genet* 2023)

---

### 2. `manuscript-writing-nature`

**Advanced scientific manuscript writing skill with journal-specific profiles, an automated 10-pass structural quality audit, and a deep iterative AI self-review engine with 12 critical lenses.**

| Metric | Value |
|---|---|
| **Lines of code/docs** | 3,711 |
| **Files** | 17 |
| **Review lenses** | 12 |
| **Exemplary papers** | 9 (Nature, Cell, Science, Nat Rev Genet, Nat Mach Intell) |
| **Output formats** | LaTeX + DOCX |

#### Pipeline Overview

```
Journal Selection → Story Arc → Figures First → Drafting → 10-Pass Audit → Deep AI Review → Formatting → Submission
```

#### What It Does

| Step | Description |
|---|---|
| **1. Journal Selection** | Scope check, 90-second test, acceptance rate assessment |
| **2. Story Arc** | Hook → Context → Gap → Approach → Finding → Impact framework |
| **3. Figures First** | Visual narrative before text, self-explanatory captions |
| **4. Section Drafting** | Methods → Results → Discussion → Intro → Abstract → Title order |
| **5a. Structural Audit** | 10-pass automated + manual checks (word count, citations, stats, figures) |
| **5b. Deep AI Review** | 12-lens iterative review until zero concerns (narrative, AI vocab, biology, clarity) |
| **6–8. Submission** | Formatting, cover letter, revision handling with rebuttal templates |

#### The 12 Review Lenses

| Lens | What It Checks |
|---|---|
| Narrative flow | Logical transitions, story momentum |
| Figure captions | Self-explanatory? Statistics included? |
| Biological interpretation | Grounded in high-impact literature? |
| AI vocabulary | Detects and eliminates 20+ AI cliché patterns |
| Clarity & precision | Unambiguous pronouns, specific numbers |
| Methods reproducibility | Could someone replicate from Methods alone? |
| Data-claim alignment | Every claim traceable to data |
| Logical consistency | Numbers match across sections |
| Impact framing | "So what?" test on every section |
| Tone & register | Confident without arrogance |
| Completeness | Abbreviations, limitations, alternatives |
| Reader experience | Would you accept as a reviewer? |

#### Reference Documentation

| Document | Coverage |
|---|---|
| `nature_guide.md` | Complete Nature author guide: structure, word limits, figure specs, references |
| `writing_style.md` | Narrative hooks, funnel structure, active voice, power words, abstract formula |
| `ai_review_protocol.md` | Deep iterative 12-lens self-review with naive reader test and convergence criteria |
| `quality_checklist.md` | 10-pass structural audit with checklists for each pass |
| `exemplary_articles.md` | 9 real top papers (AlphaFold, CMap, Scanpy, SHAP) with style lessons |
| `statistical_reporting.md` | P-values, effect sizes, ML metrics, enrichment, survival analysis |
| `figure_standards.md` | Dimensions, DPI, color-blind palettes, Extended Data vs Supplementary |
| `reviewer_expectations.md` | Top 10 rejection reasons, rebuttal format, reviewer simulation |
| `citation_formatting.md` | Nature citation style, BibTeX templates, journal abbreviations |
| `data_code_availability.md` | FAIR principles, repository requirements, template statements |

#### Automated Scripts

| Script | Purpose |
|---|---|
| `quality_check.py` | Automated quality engine — checks structure, word count, acronyms, figures |
| `scaffold_manuscript.py` | Generates complete project directory with LaTeX/DOCX templates |
| `word_counter.py` | Section-by-section word counting with journal limit comparison |
| `citation_validator.py` | Validates BibTeX against manuscript — finds missing/uncited references |

---

## 🚀 Installation

### Option 1: Install a single skill

Copy the skill folder to your Antigravity skills directory:

```bash
# Clone this repository
git clone https://github.com/Walter-Odur/Advanced_Antigravity_Skills_For_Bioinformatics_Research.git

# Copy the skill you need
cp -r Advanced_Antigravity_Skills_For_Bioinformatics_Research/scrna-drug-repurposing \
      ~/.gemini/config/skills/scrna-drug-repurposing

# Or the manuscript writing skill
cp -r Advanced_Antigravity_Skills_For_Bioinformatics_Research/manuscript-writing-nature \
      ~/.gemini/config/skills/manuscript-writing
```

**Windows (PowerShell):**
```powershell
git clone https://github.com/Walter-Odur/Advanced_Antigravity_Skills_For_Bioinformatics_Research.git

# scRNA-seq skill
Copy-Item -Recurse "Advanced_Antigravity_Skills_For_Bioinformatics_Research\scrna-drug-repurposing" `
    "$env:USERPROFILE\.gemini\config\skills\scrna-drug-repurposing"

# Manuscript writing skill
Copy-Item -Recurse "Advanced_Antigravity_Skills_For_Bioinformatics_Research\manuscript-writing-nature" `
    "$env:USERPROFILE\.gemini\config\skills\manuscript-writing"
```

### Option 2: Install all skills

```bash
git clone https://github.com/Walter-Odur/Advanced_Antigravity_Skills_For_Bioinformatics_Research.git
cp -r Advanced_Antigravity_Skills_For_Bioinformatics_Research/*/ ~/.gemini/config/skills/
```

### Python Dependencies

```bash
pip install scanpy scikit-learn xgboost lightgbm shap optuna matplotlib seaborn \
            pandas numpy scipy requests gseapy joblib
```

---

## 📖 Usage

Once installed, Antigravity will automatically activate the skill when your task matches its domain. You can also explicitly invoke it:

> *"Build a drug repurposing pipeline from my scRNA-seq dataset using the scrna-drug-repurposing skill"*

Each skill includes:
- **`SKILL.md`** — Main instruction file with step-by-step guidance
- **`references/`** — Deep-dive documentation on each pipeline component
- **`examples/`** — Complete reference implementations with real-world datasets
- **`scripts/`** — Utility scripts (e.g., project scaffolding)

---

## 🗺️ Roadmap — Upcoming Skills

| Skill | Domain | Status |
|---|---|---|
| `scrna-drug-repurposing` | scRNA-seq → ML → Drug Repurposing | ✅ **Available** |
| `manuscript-writing-nature` | Scientific manuscript writing for Nature | ✅ **Available** |
| `spatial-transcriptomics` | Spatial gene expression analysis (Visium, MERFISH, Slide-seq) | 🔜 Planned |
| `multi-omics-integration` | Integrating transcriptomics, proteomics, metabolomics | 🔜 Planned |
| `crispr-screen-analysis` | Genome-wide CRISPR screen hit calling and pathway analysis | 🔜 Planned |
| `phylogenomics` | Phylogenetic analysis, molecular evolution, tree building | 🔜 Planned |
| `protein-structure-ml` | AlphaFold-based structure prediction and drug-target analysis | 🔜 Planned |
| `metagenomics-16s` | 16S rRNA amplicon sequencing analysis pipeline | 🔜 Planned |
| `clinical-biomarker` | Biomarker discovery from clinical multi-omics data | 🔜 Planned |

---

## 📚 References

The skills in this repository are grounded in methodology from 70+ peer-reviewed publications. Key references:

### Core Pipeline
1. Subramanian A, et al. (2017). A Next Generation Connectivity Map. *Cell*, 171(6):1437–1452
2. Lundberg SM, Lee S-I. (2017). A Unified Approach to Interpreting Model Predictions. *NeurIPS 2017*
3. Lundberg SM, et al. (2020). From local explanations to global understanding. *Nat Mach Intell*, 2:56–67

### scRNA-seq
4. Heumos L, et al. (2023). Best practices for single-cell analysis across modalities. *Nat Rev Genet*, 24:550–572
5. Luecken MD, Theis FJ. (2019). Current best practices in scRNA-seq analysis. *Mol Syst Biol*, 15(6):e8746
6. Wolf F, et al. (2018). SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol*, 19:15

### Machine Learning
7. Breiman L. (2001). Random Forests. *Machine Learning*, 45:5–32
8. Chen T, Guestrin C. (2016). XGBoost. *KDD 2016*
9. Ke G, et al. (2017). LightGBM. *NeurIPS 2017*
10. Chicco D, Jurman G. (2020). MCC advantages over F1 and accuracy. *BMC Genomics*, 21:6

### Drug Repurposing
11. Pushpakom S, et al. (2019). Drug repurposing: progress and recommendations. *Nat Rev Drug Discov*, 18:41–58
12. Corsello SM, et al. (2020). Discovering anticancer potential of non-oncology drugs. *Nat Cancer*, 1:235–248
13. Lamb J, et al. (2006). The Connectivity Map. *Science*, 313(5795):1929–1935

### Manuscript Writing
14. Jumper J, et al. (2021). Highly accurate protein structure prediction with AlphaFold. *Nature*, 596:583–589
15. Gideon HP, et al. (2022). Multimodal profiling of lung granulomas in macaques. *Immunity*, 55:827–846
16. Kaufmann SHE, et al. (2018). Host-directed therapies for bacterial and viral infections. *Nat Rev Drug Discov*, 17:35–56

> See individual skill `references/` directories for the complete citation list.

---

## 👤 Author

**Walter Odur**

- GitHub: [@Walter-Odur](https://github.com/Walter-Odur)

---

*Built with [Antigravity IDE](https://github.com/google-deepmind/antigravity) — Advanced Agentic AI Coding*
