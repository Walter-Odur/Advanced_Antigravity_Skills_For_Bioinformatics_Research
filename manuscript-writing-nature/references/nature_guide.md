# Nature Journal — Complete Author Guide

## Article Format (Since 2019)

Nature retired the "Letter" format in 2019. All primary research is now published as **Articles**.

---

## Manuscript Structure

Nature Articles follow this exact order:

### 1. Title
- **Maximum**: ~90 characters (fits two lines in print)
- **Style**: Concise, informative, active voice preferred
- **Prohibited**: Abbreviations, formulae, trade names, punctuation at end
- **Examples of strong Nature titles**:
  - "Single-cell transcriptomics reveals immune signatures in tuberculosis granulomas"
  - "Machine learning identifies host-directed therapeutic targets for drug-resistant TB"
  - "A next generation connectivity map: L1000 platform and the first 1,000,000 profiles"

### 2. Authors & Affiliations
- Full names, not initials only
- Numbered affiliations
- Corresponding author(s) marked with asterisk
- ORCID IDs required for all authors
- Equal contributions noted with daggers (†)

### 3. Summary (Abstract)
- **Format**: Single, unstructured, **bolded** paragraph
- **Length**: ≤200 words
- **NO references** within the summary
- **NO section labels** (Background:, Methods:, etc.)
- **NO abbreviations** unless absolutely essential
- **Structure** (implicit, not labeled):
  1. Context — one sentence establishing broad significance
  2. Problem — one sentence identifying the gap
  3. Approach — one to two sentences on what you did
  4. Key findings — two to three sentences on main results
  5. Impact — one sentence on implications

**Template**:
```
[Broad context establishing why this matters to science]. However,
[specific knowledge gap or unresolved problem]. Here we [approach/method],
[using what tools/data]. We show that [main finding 1], and demonstrate
that [main finding 2]. [Additional finding if needed]. These results
[implication for the field], suggesting [broader impact or future direction].
```

### 4. Main Body Text

**Critical: Nature Articles typically do NOT use section headers.** The text flows as continuous paragraphs. However, the logical structure follows:

#### Introduction (~500 words)
- **Paragraph 1**: Hook + broad context (accessible to ALL scientists)
- **Paragraph 2**: Narrow to your specific field, cite key prior work (5–8 refs)
- **Paragraph 3**: Identify the gap — what remains unknown
- **Final sentence**: "Here we show/report/demonstrate that..." (your contribution)

#### Results (~1,500 words)
- Each paragraph tied to a figure panel
- Pattern: "We found X (Fig. 1a). This was associated with Y (Fig. 1b)."
- Past tense for your findings
- Present every figure panel explicitly in the text
- NO interpretation — just present findings
- Subheadings ARE allowed in Results (brief, informative)

#### Discussion (~1,000 words)
- **Paragraph 1**: Restate main finding in broader context
- **Paragraphs 2–3**: Compare to prior literature, explain discrepancies
- **Paragraph 4**: Mechanistic interpretation or proposed model
- **Paragraph 5**: Limitations (acknowledge honestly)
- **Final paragraph**: Future directions, broader implications
- Do NOT introduce new data
- Do NOT repeat Results verbatim

### 5. Methods (Online Only)
- Appears after References in print (online-only)
- **No word limit**
- Past tense throughout
- Sufficient detail for reproducibility
- Subsections with headers:
  - Study design / Data sources
  - Data processing and quality control
  - Computational methods / Models
  - Statistical analysis
  - Software and tools (with versions)
  - Ethics statement (if applicable)

### 6. References
- **Maximum**: 50 references for Articles
- **Style**: Numbered superscripts in text
- **Format**: `Author1, A. B., Author2, C. D. & Author3, E. F. Title. J. Abbrev. vol, pages (year).`
- Consecutive references: `1,2,3` or `1–3`
- Journal names abbreviated (italicized)
- All authors listed for ≤5 authors; for 6+, list first author then "et al."

### 7. End Matter (Required Statements)

Each must be a separate section:

#### Acknowledgements
- Funding sources with grant numbers
- Technical assistance
- Do NOT include authors here

#### Author Contributions
- CRediT taxonomy or free text
- Example: "W.O. conceived and designed the study, performed the computational analysis, and wrote the manuscript. X.Y. supervised the project and reviewed the manuscript."

#### Competing Interests
- Must declare or state: "The authors declare no competing interests."

#### Data Availability
- **Mandatory**
- State where data can be accessed
- Include accession numbers, DOIs, or repository URLs
- If data cannot be shared, justify why

#### Code Availability
- **Mandatory** if custom code was used
- Link to public repository (GitHub, Zenodo with DOI)
- Include version information

### 8. Extended Data
- Up to 10 figures/tables
- Peer-reviewed (unlike Supplementary Information)
- Same quality standards as main figures
- Referenced as "Extended Data Fig. 1" in text

### 9. Supplementary Information
- Not peer-reviewed to the same standard
- Large datasets, additional methods, supplementary figures/tables
- Referenced as "Supplementary Fig. 1" or "Supplementary Table 1"

---

## Word Count Rules

| Component | Limit | Counted? |
|---|---|---|
| Summary/Abstract | ≤200 words | Separate count |
| Main body text | ≤3,000 words | Yes |
| Methods | Unlimited | Not counted |
| Figure legends | ≤300 words each | Not counted |
| References | ≤50 items | Not counted |
| Acknowledgements | — | Not counted |

**What counts as "main body text"**: Introduction + Results + Discussion only. Methods, references, figure legends, and end-matter statements are excluded.

---

## Figure & Table Specifications

### Quantity
- Maximum **6 display items** (figures + tables combined) in main paper
- Maximum **10 Extended Data** items
- Unlimited Supplementary figures/tables

### Dimensions
| Type | Width |
|---|---|
| Single column | 89 mm (3.5 in) |
| 1.5 columns | 120 mm (4.7 in) |
| Double column | 183 mm (7.2 in) |
| Full page height | 247 mm (9.7 in) |

### Resolution
| Content Type | Minimum DPI |
|---|---|
| Photographs / micrographs | 300 |
| Line art / graphs | 600 |
| Combination (line art + photo) | 500 |

### Format
- **Initial submission**: JPEG, PNG, TIFF, PDF, or EPS
- **Accepted revision**: TIFF, EPS, or PDF (high-resolution)
- **Embed in manuscript** for initial submission (separate files for revision)

### Labels & Annotations
- **Font**: Arial or Helvetica
- **Panel labels**: Lowercase bold (**a**, **b**, **c**), 8–12 pt
- **Axis labels**: 8–10 pt, with units
- **Scale bars**: Include for microscopy images
- **Color-blind safe**: Avoid red-green only distinctions
  - Recommended palette: blue (#0072B2), orange (#E69F00), green (#009E73), vermillion (#D55E00), sky blue (#56B4E9), yellow (#F0E442), black (#000000)

### Figure Legends
- **Maximum**: 300 words per legend
- **Location**: End of manuscript (NOT on the figure files)
- **Structure**:
  1. **Title line** (bold): One-sentence description of what the figure shows
  2. **Panel descriptions**: Brief description of each panel (a, b, c...)
  3. **Statistical information**: Sample sizes (n), statistical tests, P values
  4. **Scale information**: Scale bars, axis descriptions
  5. **Data representation**: "Data are mean ± s.d." or "Box plots show median, IQR..."

**Example legend**:
```
**Fig. 1 | Ensemble machine learning identifies TB-associated transcriptomic signatures.**
**a**, UMAP visualization of 15,247 single cells colored by disease status (red, TB-affected;
blue, control). **b**, ROC curves for Random Forest (blue), XGBoost (green), LightGBM (orange),
and Stacking Ensemble (purple) classifiers. AUC values shown in legend. **c**, Top 20 genes
ranked by mean absolute SHAP value. Bar color indicates direction of association (red,
upregulated in TB; blue, downregulated). n = 12,198 training cells, 3,049 test cells.
Statistical significance assessed by DeLong test for AUC comparison.
```

---

## Reference Format

### In-Text Citations
- Superscript numbers: `... as shown previously¹.`
- Consecutive: `... multiple studies¹⁻³ have demonstrated...`
- Non-consecutive: `... prior work¹,⁴,⁷ suggests...`
- Position: After punctuation, not before

### Reference List Format

**Journal article** (≤5 authors — list all):
```
1. Subramanian, A., Narayan, R., Corsello, S. M., Peck, D. D. & Bhatt, D. L.
   A next generation connectivity map: L1000 platform and the first 1,000,000
   profiles. Cell 171, 1437–1452 (2017).
```

**Journal article** (6+ authors — first author et al.):
```
2. Heumos, L. et al. Best practices for single-cell analysis across modalities.
   Nat. Rev. Genet. 24, 550–572 (2023).
```

**Book**:
```
3. Hastie, T., Tibshirani, R. & Friedman, J. The Elements of Statistical
   Learning 2nd edn (Springer, 2009).
```

**Book chapter**:
```
4. Smith, J. in Computational Biology (ed. Jones, K.) 45–78 (Academic Press, 2020).
```

**Preprint**:
```
5. Doe, J. et al. Novel approach to drug repurposing. Preprint at
   https://doi.org/10.1101/2024.01.15.12345 (2024).
```

### Journal Abbreviations
Use standard MEDLINE/PubMed abbreviations:
- *Nature* → `Nature`
- *Nature Methods* → `Nat. Methods`
- *Nature Machine Intelligence* → `Nat. Mach. Intell.`
- *Cell* → `Cell`
- *Science* → `Science`
- *Genome Biology* → `Genome Biol.`
- *BMC Genomics* → `BMC Genomics`
- *Bioinformatics* → `Bioinformatics`

---

## Required Reporting Summary

Nature requires a completed **Reporting Summary** form:
- Download from: https://www.nature.com/documents/nr-reporting-summary.pdf
- Covers: statistics, data collection, study design, code availability
- Must be submitted alongside the manuscript

---

## Submission Format

### Initial Submission
- **File format**: PDF, Word (.docx), or LaTeX (.tex + .bib)
- **Figures**: Embedded in manuscript OR separate files
- **Line numbers**: Required
- **Double spacing**: Required
- **Font**: 12-pt standard font (Times New Roman, Arial)

### Revision (After Peer Review)
- **File format**: Word with track changes OR LaTeX with diff
- **Figures**: Separate high-resolution files
- **Response to reviewers**: Point-by-point rebuttal document
- **Reporting Summary**: Updated if needed

---

## Cover Letter Template

```
Dear Dr. [Editor Name],

We submit our manuscript entitled "[Title]" for consideration as an
Article in Nature.

[1-2 sentences summarizing the main finding and its significance,
written for a non-specialist.]

[1 sentence explaining why this work is appropriate for Nature's
broad readership.]

[1 sentence on the methodology used.]

This work has not been published elsewhere and is not under
consideration by another journal. All authors have approved the
manuscript and agree with its submission to Nature.

We suggest the following reviewers:
1. [Name, Institution, Email] — expert in [relevant area]
2. [Name, Institution, Email] — expert in [relevant area]
3. [Name, Institution, Email] — expert in [relevant area]

We would prefer to exclude [Name] due to [brief reason].

Sincerely,
[Corresponding Author Name]
[Institution]
[Email]
[ORCID]
```
