# Nature Narrative Writing Style Guide

## The Nature Voice

Nature papers are fundamentally different from papers in specialist journals. Nature demands **accessible, story-driven prose** that captivates scientists across disciplines. A cell biologist should find your machine learning paper compelling. A physicist should understand why your genomics discovery matters.

---

## The 90-Second Test

Nature editors triage manuscripts in under 90 seconds. In that time, they read:
1. The **title** (~5 seconds)
2. The **summary/abstract** (~30 seconds)
3. **Figure 1** and its legend (~30 seconds)
4. The **first paragraph** of the main text (~25 seconds)

If the significance is not crystal clear from these four elements, the manuscript is desk-rejected.

**Self-test**: Cover everything except these four elements. Can a non-specialist grasp the core advance? If not, rewrite them.

---

## The Hook: First Sentence

The opening sentence is the most important sentence in the paper. It must:

1. **Capture a broad audience** — no jargon, no field-specific assumptions
2. **Establish significance** — why should anyone care?
3. **Create tension** — set up a problem that needs solving

### Good Hooks (Nature-style)

| Hook Type | Example |
|---|---|
| **Significance hook** | "Tuberculosis remains the leading infectious disease killer worldwide, claiming 1.3 million lives annually." |
| **Paradox hook** | "Despite decades of antibiotic development, drug-resistant tuberculosis continues to outpace therapeutic innovation." |
| **Discovery hook** | "Single-cell technologies have revealed unexpected heterogeneity in immune responses to infection." |
| **Numbers hook** | "More than 10 million people develop active tuberculosis each year, yet fewer than 5% receive adequate treatment." |

### Bad Hooks (Avoid)

| Problem | Example | Why It Fails |
|---|---|---|
| Too narrow | "The role of CCL18 in granuloma formation..." | Only specialists care |
| Too technical | "Stacking ensemble classifiers with L1-penalized..." | Jargon in first sentence |
| Too passive | "It has been shown that tuberculosis..." | Weak, indirect |
| Too vague | "Machine learning is increasingly used in biology." | Generic, no tension |

---

## The Funnel Structure (Introduction)

Nature introductions follow a **funnel** — broad to narrow:

```
┌─────────────────────────────────┐  ← Broad: significance (all scientists)
│  Why does this problem matter?  │
├───────────────────────────┤       ← Narrower: field context
│  What has been done?      │
├─────────────────────┤             ← Narrow: specific gap
│  What's missing?    │
├───────────────┤                   ← Precise: your contribution
│  Here we...   │
└───────────────┘
```

### Paragraph-by-Paragraph Blueprint

**Paragraph 1 — The World (5–6 sentences)**:
- Sentence 1: Hook (broad significance, no jargon)
- Sentence 2–3: Context establishing why this problem matters
- Sentence 4–5: Current understanding at a high level
- Sentence 6: Transition to the specific field

**Paragraph 2 — The Field (4–5 sentences)**:
- Prior work in your specific area (cite 5–8 key papers)
- What has been achieved so far
- What tools/approaches exist
- Logical progression toward the gap

**Paragraph 3 — The Gap (3–4 sentences)**:
- "However, ..." or "Despite these advances, ..."
- What remains unknown or unresolved
- Why existing approaches are insufficient
- The specific question your study addresses

**Paragraph 4 — Your Contribution (2–3 sentences)**:
- "Here we show/report/demonstrate that..."
- Brief overview of approach
- Preview of main finding (1 sentence)

### Transition Words That Work

| Purpose | Words |
|---|---|
| Building on prior work | "Building on this", "Extending these findings" |
| Contrasting | "However", "In contrast", "Despite this" |
| Advancing the argument | "Moreover", "Furthermore", "In addition" |
| Consequence | "As a result", "Consequently", "These findings suggest" |
| Temporal sequence | "Subsequently", "Following this", "More recently" |

---

## Results Writing

### The Figure-First Principle

Every Results paragraph should be anchored to a figure panel. The pattern:

```
[Context/setup sentence]. We [action] and found that [finding]
(Fig. Xa). [Elaboration on the finding]. [Quantitative detail].
This [result] was consistent with [comparison/expectation]
(Fig. Xb). [Additional observation].
```

### Active vs. Passive Voice

Nature strongly prefers **active voice**:

| ❌ Passive (Avoid) | ✅ Active (Preferred) |
|---|---|
| "It was found that CCL18 was upregulated" | "We found that CCL18 was upregulated" |
| "The analysis was performed using..." | "We performed the analysis using..." |
| "A significant difference was observed" | "We observed a significant difference" |
| "The model was trained on..." | "We trained the model on..." |

### Tense Rules

| Section | Tense | Example |
|---|---|---|
| Introduction (established facts) | Present | "TB **is** the leading cause..." |
| Introduction (prior studies) | Past | "Smith et al. **showed** that..." |
| Results (your findings) | Past | "We **found** that..." |
| Results (referring to figures) | Present | "Fig. 1a **shows** the distribution..." |
| Discussion (your results) | Past | "Our analysis **revealed**..." |
| Discussion (implications) | Present | "These findings **suggest**..." |
| Methods | Past | "We **used** Scanpy v1.9..." |

---

## Discussion Architecture

### The Inverted Funnel

The Discussion mirrors the Introduction — but inverted (narrow to broad):

```
└───────────────┘                   ← Narrow: your main finding
├─────────────────────┤             ← Broader: comparison to literature
├───────────────────────────┤       ← Broader: mechanistic model
├─────────────────────────────────┤ ← Broadest: implications + future
┌─────────────────────────────────┐
```

### Paragraph-by-Paragraph Blueprint

**Paragraph 1 — Restate the Advance (3–4 sentences)**:
- "In this study, we demonstrate that..."
- Main finding in broader context
- How this changes our understanding

**Paragraphs 2–3 — Literature Comparison (4–5 sentences each)**:
- How your findings relate to prior work
- Agreements: "Consistent with X, we found..."
- Disagreements: "In contrast to Y, our results suggest..."
- Explain discrepancies (different methods, populations, etc.)

**Paragraph 4 — Mechanistic Interpretation (3–4 sentences)**:
- Proposed mechanism or model
- "One possible explanation is..."
- Link molecular findings to biological function
- This is where you can speculate (carefully)

**Paragraph 5 — Limitations (3–4 sentences)**:
- Acknowledge honestly — reviewers WILL find them
- Sample size, generalizability, technical limitations
- "While our study provides X, several limitations should be noted..."
- Frame limitations constructively: "Future studies with larger cohorts could..."

**Paragraph 6 — Future Directions & Impact (3–4 sentences)**:
- What should be done next
- Broader implications for the field
- Practical/clinical applications
- End on a strong, forward-looking note

---

## Sentence-Level Craft

### Sentence Length

- **Target**: 15–25 words per sentence on average
- **Vary length**: Mix short punchy sentences with longer complex ones
- **Flag**: Any sentence >40 words should be split
- **Rule**: One idea per sentence

### Paragraph Length

- **Target**: 4–6 sentences per paragraph
- **Maximum**: 8 sentences (split longer paragraphs)
- **Minimum**: 2 sentences (single-sentence paragraphs feel fragmented)

### Power Words for Nature Papers

| Instead of... | Use... |
|---|---|
| "We looked at" | "We investigated / examined / characterized" |
| "We used" | "We employed / applied / utilized" |
| "Big" | "Substantial / considerable / marked" |
| "Important" | "Critical / essential / pivotal" |
| "Interesting" | "Noteworthy / striking / remarkable" |
| "Shows" | "Reveals / demonstrates / establishes" |
| "Suggests" | "Indicates / implies / points to" |
| "A lot of" | "Numerous / substantial / extensive" |

### Words to Avoid

| Word | Why | Better Alternative |
|---|---|---|
| "Very" | Weak intensifier | Remove or use specific qualifier |
| "Basically" | Filler | Delete |
| "Obviously" | Condescending | "Notably" or delete |
| "Novel" | Overused, let the work speak | "Previously uncharacterized" |
| "Firstly" | Unnecessary formal | "First" |
| "In order to" | Wordy | "To" |
| "Due to the fact that" | Wordy | "Because" |
| "It is worth noting that" | Throat-clearing | Delete, just state the point |

---

## Common Mistakes That Get Papers Rejected

### 1. Starting Too Narrow
❌ "CCL18, a chemokine produced by alternatively activated macrophages, has been implicated in fibrotic lung diseases..."
✅ "The host immune response to tuberculosis involves complex cellular interactions that determine clinical outcome..."

### 2. Results That Read Like a Lab Notebook
❌ "We then ran XGBoost with max_depth=6 and learning_rate=0.1. The AUC was 0.943."
✅ "Gradient-boosted ensemble methods achieved robust classification of TB-affected cells (AUC = 0.943; Fig. 2b), outperforming individual classifiers."

### 3. Discussion That Repeats Results
❌ "Our results showed that CCL18 was the most important gene (SHAP = 0.023). S100A12 was the second most important gene..."
✅ "The identification of CCL18 as the dominant classifier of TB disease status aligns with its established role in alternatively activated macrophage biology and suggests a previously unrecognized function in granuloma-associated immune regulation."

### 4. Passive Throughout
❌ "The data were analyzed and it was found that a significant difference was present between groups."
✅ "We analyzed the single-cell data and found a significant difference between TB-affected and control cells."

### 5. Missing the "So What?"
Every finding needs interpretation. After presenting a result, answer:
- **So what?** — Why does this matter?
- **Compared to what?** — How does it relate to prior work?
- **What next?** — What question does this raise?

---

## The Abstract Formula (200 Words)

Allocate words approximately:

| Component | Words | Content |
|---|---|---|
| Context | ~30 | Broad significance + specific problem |
| Gap | ~20 | What remains unknown |
| Approach | ~30 | What you did (methods, briefly) |
| Findings | ~80 | Main results (quantitative) |
| Impact | ~40 | Implications + significance |
| **Total** | **~200** | |

### Example (TB Drug Repurposing):

> Tuberculosis remains a leading cause of death from a single infectious agent, and the emergence of drug-resistant strains demands new therapeutic strategies. Host-directed therapies, which modulate the host immune response rather than targeting the pathogen directly, represent a promising but underexplored approach. Here we apply an ensemble machine learning framework—integrating Random Forest, XGBoost, LightGBM, and a stacking classifier—to single-cell RNA sequencing data from 15,247 cells across tuberculosis-affected and control lung tissue. We identify a 127-gene disease signature through interventional SHAP analysis, with CCL18, S100A12, and COL3A1 emerging as the most discriminative features. Using this signature for connectivity map-based drug repurposing, we identify 23 compounds that significantly reverse the TB transcriptomic profile, including several FDA-approved drugs with established safety profiles. Our integrative approach bridges single-cell genomics, interpretable machine learning, and computational pharmacology to nominate host-directed therapeutic candidates for tuberculosis.

(198 words)

---

## DOCX-Specific Formatting Notes

When outputting manuscripts in DOCX format:

### Style Requirements
- **Font**: 12-pt Times New Roman or Arial
- **Line spacing**: Double-spaced throughout
- **Margins**: 1 inch (2.54 cm) all sides
- **Page numbers**: Bottom center or bottom right
- **Line numbers**: Continuous throughout (use Word's line numbering feature)

### Section Formatting
- Title: 14-pt bold, centered
- Author names: 12-pt, centered
- Summary: Bold paragraph, left-aligned
- Main text: 12-pt, left-aligned, with paragraph indentation (0.5 in)
- References: Numbered list, hanging indent
- Figure legends: Grouped at end, each starting with bold title

### Track Changes for Revision
- Use Word's built-in track changes for revisions
- Accept all changes in the clean version
- Provide both tracked and clean copies
