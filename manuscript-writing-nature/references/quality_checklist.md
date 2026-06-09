# 10-Pass Iterative Quality Audit Protocol

## Overview

This document defines the structured 10-pass quality audit system for scientific manuscripts. Each pass has a single focus area. **All issues found in a pass must be fixed before moving to the next pass.** After completing all 10 passes, re-run any pass that required significant changes.

The philosophy: **separate concerns, fix one thing at a time.** Trying to check everything at once leads to missed errors.

---

## Pass 1: Structure & Completeness

**Focus**: Are all required sections present in the correct order?

### Checklist

```
□ Title present and ≤90 characters
□ Author list with affiliations
□ Summary/Abstract present
□ Main body text present
□ Methods section present
□ References section present
□ Acknowledgements present
□ Author Contributions statement
□ Competing Interests declaration
□ Data Availability statement
□ Code Availability statement (if applicable)
□ Extended Data referenced (if applicable)
□ Supplementary Information referenced (if applicable)
□ Sections in correct journal order
□ No orphaned sections or placeholders
□ No TODO/FIXME comments remaining
```

### Common Issues
- Missing Data Availability statement (automatic desk rejection at Nature)
- Missing Code Availability (required if custom code was used)
- Author Contributions missing or incomplete
- Sections in wrong order (Methods before Results, etc.)

---

## Pass 2: Word Count Compliance

**Focus**: Does every section meet the journal's word limits?

### Nature Limits

| Section | Limit | Action if Over |
|---|---|---|
| Summary/Abstract | ≤200 words | Tighten language, remove qualifiers |
| Main body (Intro + Results + Discussion) | ≤3,000 words | Move detail to Methods or Supplementary |
| Figure legends | ≤300 words each | Condense statistical details |
| Title | ≤90 characters | Shorten, remove unnecessary words |

### Automated Check
```bash
python scripts/word_counter.py manuscript/main.tex --profile nature
```

### Reduction Strategies
1. **Remove throat-clearing**: "It is worth noting that..." → delete
2. **Remove redundancies**: "completely eliminated" → "eliminated"
3. **Compress citations**: "Smith (2020) showed that X. Jones (2021) confirmed X." → "X has been confirmed¹²."
4. **Move to Methods**: Detailed experimental descriptions belong in Methods
5. **Move to Supplementary**: Secondary analyses, validation experiments

---

## Pass 3: Narrative Flow & Story Arc

**Focus**: Does the manuscript tell a compelling, logical story?

### Checklist

```
□ First sentence hooks a broad audience
□ Introduction follows funnel structure (broad → narrow)
□ Gap is clearly identified before the contribution statement
□ "Here we..." sentence clearly states the contribution
□ Results follow a logical sequence (not just chronological)
□ Each Results paragraph is anchored to a figure
□ Discussion opens by restating the main finding
□ Discussion does NOT repeat Results verbatim
□ Discussion includes literature comparison
□ Limitations are honestly acknowledged
□ Final paragraph looks forward (implications, future work)
□ Red thread: can you follow the argument from title → abstract → conclusion?
```

### The "Paragraph Topic Sentence" Test
Read ONLY the first sentence of every paragraph. They should form a logical argument on their own:

```
1. "Tuberculosis remains a leading cause of death..." (hook)
2. "Single-cell RNA sequencing has emerged..." (context)
3. "However, computational approaches for..." (gap)
4. "Here we apply an ensemble ML framework..." (contribution)
5. "We first quality-controlled and preprocessed..." (Results start)
6. "Ensemble classifiers achieved high accuracy..." (key finding)
7. "SHAP analysis identified 127 genes..." (signature)
8. "Drug repurposing identified 23 candidates..." (application)
9. "Our integrative approach demonstrates..." (Discussion start)
10. "The identification of CCL18 as..." (interpretation)
11. "Several limitations should be noted..." (limitations)
12. "These findings establish a framework..." (outlook)
```

If the topic sentences don't form a coherent argument, restructure.

---

## Pass 4: Accessibility & Jargon

**Focus**: Can a scientist outside your subfield understand the paper?

### Checklist

```
□ First paragraph contains NO jargon or field-specific acronyms
□ All acronyms defined at first use
□ Technical terms explained or contextualized
□ No undefined abbreviations
□ Figures understandable without reading the text
□ Summary/Abstract uses no acronyms (or minimal, defined ones)
□ Methods clearly describes tools for non-specialists
□ No insider references ("as is well known..." — to whom?)
```

### Jargon Red Flags
Test: Could a PhD student in a different field understand each sentence?

| Too Jargon-Heavy | Accessible Alternative |
|---|---|
| "We computed SHAP values using TreeExplainer with interventional feature perturbation" | "We used SHAP (SHapley Additive exPlanations), an interpretable machine learning framework, to quantify each gene's contribution to classification" |
| "The StackingClassifier used L1-penalized logistic regression as meta-learner" | "We combined predictions from multiple classifiers using a stacking ensemble approach" |
| "HVG selection via Seurat v3 method" | "We selected the 2,000 most variable genes across cells" |

---

## Pass 5: Figures & Tables

**Focus**: Are display items compliant, clear, and self-explanatory?

### Checklist

```
□ Total display items ≤6 (figures + tables)
□ Every figure referenced in text ("Fig. 1a")
□ Every table referenced in text ("Table 1")
□ Figure legends ≤300 words each
□ Legends include: title, panel descriptions, sample sizes, statistics
□ Panel labels are lowercase bold (a, b, c)
□ Font is Arial/Helvetica, 8-12 pt
□ Colors are color-blind safe
□ Resolution ≥300 DPI (≥600 for line art)
□ Axis labels include units
□ Scale bars included for microscopy
□ No redundancy between figures and tables
□ Figures tell the story independently of text
□ Extended Data items ≤10 (if used)
```

### Figure Quality Test
Print the figures at journal column width (89mm or 183mm). Can you read all labels? If not, the font is too small.

---

## Pass 6: Statistical Reporting

**Focus**: Are all statistical claims properly supported?

### Checklist

```
□ Every p-value has an associated test name
□ Sample sizes (n) reported for all analyses
□ Effect sizes reported alongside p-values
□ Multiple testing corrections applied (Bonferroni, FDR) where needed
□ Confidence intervals reported for key estimates
□ Exact p-values given (not just "p < 0.05") unless very small
□ For very small p-values: "P < 1 × 10⁻¹⁵" not "P = 0.000"
□ Error bars defined (s.d., s.e.m., 95% CI)
□ Software and versions for all statistical analyses
□ No "trends toward significance" (p = 0.07 is not significant)
□ Cross-validation results include mean ± s.d.
□ AUC values include confidence intervals
□ All statistical methods described in Methods section
```

### Nature Statistical Requirements
- Report **exact P values** to two significant figures for P ≥ 0.001
- Report as "P < 0.001" for smaller values, or exact (e.g., P = 3.2 × 10⁻⁸)
- Always state the **statistical test** used
- Always state **n** (biological replicates, not technical replicates)
- Define all error bars in figure legends
- Report effect sizes for primary outcomes

### Common Statistical Errors
1. Using t-tests on non-normal data → use Mann-Whitney U
2. Multiple comparisons without correction → use Bonferroni or BH-FDR
3. Confusing s.d. with s.e.m. → use s.d. for describing variability, s.e.m. for comparing means
4. Reporting "significant" without p-value → always include exact values
5. n too small for the test → power analysis or note as limitation

---

## Pass 7: Citation Integrity

**Focus**: Are all citations accurate, complete, and properly formatted?

### Checklist

```
□ All factual claims have citations
□ No "citation needed" gaps
□ References ≤50 for Nature Articles
□ Mix of recent work (last 5 years) and foundational papers
□ No excessive self-citation (>20% is a flag)
□ Citations are primary sources, not secondary references
□ All references formatted in Nature style (numbered superscripts)
□ Journal names properly abbreviated
□ No broken DOIs or URLs
□ All referenced papers actually support the claim made
□ Review articles cited only for broad context, not specific findings
□ No predatory journal citations
□ Preprints noted as "Preprint at..." with DOI
```

### Automated Check
```bash
python scripts/citation_validator.py manuscript/references.bib manuscript/main.tex
```

---

## Pass 8: Required Statements & Compliance

**Focus**: Are all mandatory sections and statements present?

### Checklist

```
□ Data Availability: specific repository, accession numbers, DOIs
□ Code Availability: GitHub/Zenodo link, version, DOI
□ Author Contributions: all authors listed with roles
□ Competing Interests: declared or "none"
□ Acknowledgements: funding sources with grant numbers
□ Ethics: IRB approval if human data, IACUC if animal data
□ Informed Consent: if human subjects
□ Nature Reporting Summary: completed (downloadable form)
□ ORCID IDs: listed for all authors
□ Corresponding author: email provided
```

### Data Availability Examples

**Good**:
> The single-cell RNA-seq data used in this study are available from the Gene Expression Omnibus under accession number GSE178837. Processed data and intermediate analysis files are available at Zenodo (https://doi.org/10.5281/zenodo.XXXXXXX). Source data for all figures are provided with this paper.

**Bad**:
> Data are available upon reasonable request.
(Nature no longer accepts this for most data types)

---

## Pass 9: Consistency

**Focus**: Is the manuscript internally consistent?

### Checklist

```
□ Same term for same concept throughout (no synonym switching)
□ Gene names: consistent italics (*CCL18*, not CCL18 and CCL18)
□ Protein names: consistent non-italic (p53, not *p53*)
□ Species names: italicized (*Mycobacterium tuberculosis*)
□ Abbreviations: defined at first use, used consistently after
□ Numbers: consistent style (words for <10, digits for ≥10)
□ Tense: consistent within each section
□ Figure references: consistent format ("Fig. 1a" not "Figure 1A")
□ P-value format: consistent ("P = 0.023" or "p = 0.023", not both)
□ Decimal places: consistent precision across similar values
□ Colour/color spelling: consistent (American English for Nature)
□ Table formatting: consistent across all tables
□ Reference format: consistent throughout
□ Hyphenation: consistent (single-cell or single cell, not both)
```

### American vs. British English
Nature accepts either, but must be **consistent throughout**:
- Pick one and stick with it
- Common traps: "colour" vs. "color", "analyse" vs. "analyze", "behaviour" vs. "behavior"

---

## Pass 10: Reviewer Simulation

**Focus**: Anticipate and address objections before submission.

### Think Like a Hostile Reviewer

For each major claim, ask:

```
1. "What's the evidence for this claim?"
   → Is it supported by data in the paper?

2. "Could this result be an artifact?"
   → Are there confounders, batch effects, technical artifacts?

3. "What alternative explanations exist?"
   → Are they addressed in the Discussion?

4. "Is this generalizable?"
   → Does the paper acknowledge scope limitations?

5. "Why should I believe the statistics?"
   → Are the methods appropriate for the data type?

6. "What's missing?"
   → Are there obvious experiments or analyses NOT done?

7. "Is this truly novel?"
   → Is the advance clearly stated relative to prior work?
```

### Common Reviewer Objections & Preemptive Responses

| Objection | Where to Address |
|---|---|
| "Sample size is too small" | Methods (justify) + Limitations |
| "No validation cohort" | Limitations + Future Directions |
| "Overfitting concern" | Methods (cross-validation) + Results (held-out test) |
| "Correlative, not causal" | Discussion (acknowledge, propose validation) |
| "Only one dataset" | Limitations + suggest replication |
| "Statistical method inappropriate" | Methods (justify choice with citation) |
| "Missing comparison to state-of-art" | Results (include baseline comparisons) |
| "Clinical relevance unclear" | Discussion (translational implications) |

### The "Weakness Inoculation" Strategy

Address weaknesses before reviewers raise them:
- "While our study is limited to a single dataset, the use of stratified cross-validation and held-out testing..."
- "Although in vitro validation is beyond the scope of this computational study, we note that..."
- "We acknowledge that single-cell dissociation may alter gene expression profiles; however..."

This shows maturity and awareness, and often prevents the weakness from becoming a rejection criterion.

---

## Post-Audit Checklist

After completing all 10 structural passes:

```
□ All Pass 1-10 issues fixed
□ Re-run any pass that required major changes
□ MANDATORY: Complete the Deep Iterative AI Self-Review (Step 5b)
  □ Read entire manuscript through all 12 review lenses
  □ Minimum 3 complete iterations performed
  □ Naive reader test passed
  □ No AI vocabulary patterns remaining
  □ All figure captions self-explanatory
  □ All biological claims grounded in literature
  □ Final iteration produced ZERO new concerns
□ LaTeX compiles cleanly with no warnings
□ PDF renders correctly (figures, references, formatting)
□ Final word count within limits
□ Manuscript read aloud for flow
□ Given to colleague for 90-second test
□ Cover letter drafted
□ Submission checklist complete
□ Ready for submission
```

