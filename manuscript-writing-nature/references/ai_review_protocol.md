# Deep Iterative AI Self-Review Protocol

## Overview

This document defines the **deep iterative self-review protocol** — the most critical quality assurance step in manuscript preparation. After the structural 10-pass audit is complete, the AI performs **unlimited iterative deep reviews** until no further concerns remain.

This is NOT a checklist. This is a **thinking process** — a systematic, critical reading of the entire manuscript through multiple lenses, each time asking: "If I were reading this for the first time with zero context, would I understand everything clearly?"

---

## The Core Principle: The Naive Reader Test

> **"Critically and deeply read through the final manuscript. Assuming you are an absolutely new reader with no idea how the analysis was performed — would you understand easily, without any difficulty?"**

Every sentence, every figure, every claim must pass this test. If any part requires insider knowledge, specialized jargon, or prior context that isn't provided in the manuscript itself, it MUST be rewritten.

---

## Iteration Protocol

```
REPEAT {
    1. Read the ENTIRE manuscript from title to references
    2. Apply ALL review lenses (see below)
    3. Document every concern found
    4. Fix ALL concerns
    5. Re-read the fixed sections in context
} UNTIL no new concerns are found in a complete read-through
```

**Minimum iterations**: 3 (even if no issues are found, read 3 times)
**Maximum iterations**: No limit — continue until completely clean
**Exit condition**: A full read-through produces ZERO new concerns

### When to Request User Feedback

Pause the iteration loop and present findings to the user when:
- **Major structural changes** are needed (reordering sections, adding/removing figures)
- **Biological interpretation questions** arise that require domain expertise
- **Ambiguous claims** where the intended meaning is unclear
- **After every 3 iterations** — show the user what was fixed and what remains
- **Conflicting guidance** — when two improvements seem to conflict

Present the user with:
1. A summary of issues found in the current iteration
2. The specific fixes proposed (with before/after examples)
3. Any questions that require their domain knowledge
4. An honest assessment of current manuscript quality (scale: draft → good → excellent → submission-ready)

---

## Review Lenses

Each iteration applies these 12 review lenses simultaneously. The reviewer must actively think through each lens while reading.

### Lens 1: Narrative Flow & Coherence

Read the manuscript as a story. Does it flow like a compelling narrative?

**Questions to ask at every paragraph transition:**
- Does this paragraph logically follow from the previous one?
- Is there a clear transition connecting these ideas?
- Could the reader predict what comes next?
- Does the argument build momentum toward the conclusion?
- Is there any point where the reader would feel "lost" or "confused"?

**Red flags:**
- Abrupt topic changes without transition sentences
- Circular arguments (saying the same thing in different words)
- Missing logical steps (conclusion doesn't follow from evidence)
- "Orphan paragraphs" that don't connect to anything before or after
- The Discussion merely restates Results instead of interpreting them

**Fix pattern:**
```
BEFORE: "CCL18 was identified as the top gene. We next performed drug repurposing."
AFTER:  "The identification of CCL18 as the most discriminative gene — a known marker
         of alternatively activated macrophages — suggested that the TB disease
         signature captures meaningful immunological signals. Building on this
         validated signature, we next queried the Connectivity Map to identify
         compounds that could reverse the TB-associated transcriptomic profile."
```

### Lens 2: Figure Captions — Are They Self-Explanatory?

Every figure caption must tell a **complete mini-story** that a reader can understand WITHOUT reading the main text.

**The caption must answer:**
1. **What is shown?** — "UMAP visualization of 15,247 cells..."
2. **How was it generated?** — "...colored by disease status after Harmony integration..."
3. **What should the reader see?** — "...showing clear separation between TB-affected (red) and control (blue) cells..."
4. **What are the statistics?** — "n = 12,198 training, 3,049 test cells. AUC by DeLong test."
5. **What do annotations mean?** — "Error bars represent mean ± s.d. across 5-fold CV."

**Bad caption:**
```
Figure 1. UMAP plot of cells.
```

**Good caption:**
```
Figure 1 | Single-cell landscape reveals distinct transcriptomic states in
TB-affected lung tissue. a, UMAP visualization of 15,247 cells from TB-affected
(n = 5 donors) and control (n = 5 donors) lung tissue, colored by disease status.
Cells were integrated using Harmony to remove batch effects, and the top 2,000
highly variable genes were used for dimensionality reduction. b, Same UMAP
colored by cell type annotation (Leiden clustering, resolution = 0.8; annotated
using known marker genes; see Methods). c, Proportion of each cell type in
TB-affected versus control tissue. Error bars represent 95% CI from bootstrap
resampling (n = 1,000 iterations). *P < 0.05, **P < 0.01, ***P < 0.001,
Fisher's exact test with Benjamini-Hochberg correction.
```

**Check every figure for:**
```
□ Title sentence (bold) clearly states what the figure demonstrates
□ Every panel (a, b, c) is individually described
□ Sample sizes (n) are stated
□ Statistical tests are named
□ P values or significance indicators are included
□ Error bars/uncertainty measures are defined
□ Color coding is explained
□ Axis labels have units
□ The caption can be understood without reading the main text
□ Scale bars are included for microscopy/imaging
```

### Lens 3: Biological Interpretation — Is It Grounded in Literature?

Every biological claim must be:
1. **Supported by your data** — point to the specific figure/statistic
2. **Contextualized with literature** — cite high-impact papers that support or contrast
3. **Biologically plausible** — does the interpretation make mechanistic sense?

**Questions to ask for every biological claim:**
- Is this interpretation the most parsimonious explanation of the data?
- Are there alternative biological explanations I haven't considered?
- Would a domain expert agree with this interpretation?
- Am I over-interpreting correlative data as causal?
- Does this finding make sense in the context of known biology?
- Have I cited the most authoritative (high-impact) papers on this topic?

**Citation hierarchy for biological claims:**
1. **Primary evidence**: Cell, Nature, Science papers with direct experimental evidence
2. **Review articles**: Nature Reviews, Annual Reviews for context
3. **Foundational work**: Seminal papers that established the field
4. **Recent findings**: Last 3 years for current state of knowledge
5. **Your own data**: The figures and statistics in this paper

**Example of strong biological interpretation:**
```
WEAK: "CCL18 was the most important gene in our classifier."

STRONG: "CCL18 emerged as the most discriminative feature (mean |SHAP| = 0.023),
consistent with its established role as a chemokine secreted by alternatively
activated macrophages in chronic inflammatory environments¹². In tuberculosis,
CCL18-producing macrophages have been shown to accumulate in granulomas³ and
correlate with disease severity⁴, suggesting that our computational signature
captures a biologically meaningful axis of host immune polarization. Notably,
CCL18 has been proposed as a therapeutic target in fibrotic diseases⁵, raising
the possibility that modulating this pathway could benefit TB patients."
```

### Lens 4: AI Vocabulary Detection

AI-generated text has characteristic patterns. **Actively scan for and eliminate these:**

**AI vocabulary red flags — MUST be removed or replaced:**

| AI Pattern | Why It's Bad | Better Alternative |
|---|---|---|
| "delve into" | AI cliché | "investigate", "examine", "explore" |
| "leveraging" (overused) | AI cliché | "using", "applying", "employing" |
| "in the realm of" | Unnecessarily flowery | "in" |
| "it is worth noting that" | Throat-clearing | Delete — just state the point |
| "a testament to" | AI cliché | Rewrite the sentence |
| "underscores the importance" | Overused | "highlights", "demonstrates" |
| "a nuanced understanding" | AI cliché | Be specific about what understanding |
| "sheds light on" | AI cliché | "reveals", "demonstrates" |
| "a paradigm shift" | Hyperbolic | Be specific about what changed |
| "a multifaceted approach" | Vague | Describe the actual approach |
| "tapestry of" | AI cliché | Delete entirely |
| "a pivotal role" | Overused | "a key role" or be specific |
| "the landscape of" | AI cliché | Delete or be specific |
| "in conclusion" | Unnecessary | The reader knows it's the end |
| "it's important to note" | Filler | Delete — just state it |
| "furthermore" (every paragraph) | Monotonous | Vary transitions |
| "robust" (overused) | AI cliché when excessive | "reliable", "consistent", or be specific |
| "comprehensive" (self-praise) | Let the work speak | Delete or be specific |
| "interplay between" | AI cliché | "interaction between", "relationship between" |
| "serves as a" | Wordy | "is a" |

**The sniff test:** Read every sentence and ask: "Does this sound like something a human scientist would naturally write, or does it sound like it was generated by an AI?" If the latter, rewrite it.

**Additional AI patterns to catch:**
- Lists of three adjectives ("innovative, comprehensive, and robust")
- Excessive hedging ("It could potentially be argued that...")
- Overly formal phrasing where simple language works
- Uniform sentence structure (every sentence starts the same way)
- Missing specificity — vague claims that sound impressive but say nothing

### Lens 5: Clarity & Precision

Every sentence must communicate exactly one idea with maximum clarity.

**Questions for every sentence:**
- Can this sentence be misinterpreted?
- Is every pronoun unambiguous? ("It" — what does "it" refer to?)
- Are there unnecessary words? Can this be shorter without losing meaning?
- Is the sentence structure clear? (Subject → verb → object)
- Would a reader need to re-read this sentence? If so, rewrite it.

**Precision checklist:**
```
□ No vague pronouns ("this", "it", "these" without clear antecedent)
□ No dangling modifiers
□ No run-on sentences (>40 words — split them)
□ No double negatives
□ Numbers are specific (not "a large number" but "15,247 cells")
□ Comparisons are complete ("X was higher than Y" — by how much?)
□ Temporal claims are specific ("recently" → "in 2023")
```

### Lens 6: Methods Reproducibility

Read the Methods section and ask: **"Could I reproduce this analysis from scratch using ONLY the information in this section?"**

**For computational studies, verify:**
```
□ Input data: exact source, accession number, format
□ Software: name, version, URL
□ Parameters: every non-default parameter listed
□ Random seeds: specified for reproducibility
□ Hardware: if relevant (GPU type, memory)
□ Preprocessing: every step described in order
□ Statistical tests: named with justification
□ Thresholds: every cutoff stated (FDR < 0.05, |log2FC| > 1, etc.)
□ Model training: train/test split, cross-validation scheme, metrics
□ Feature selection: criteria and number of features
□ Visualization: tools and parameters used
```

### Lens 7: Data-Claim Alignment

For every claim in the manuscript, verify that the supporting data actually supports it.

**The tracing exercise:**
```
For each claim:
1. Find the claim in the text
2. Find the figure/table/statistic that supports it
3. Verify the figure/table actually shows what the text says
4. Verify the statistics are correctly reported
5. Verify the interpretation matches the data (not over-claimed)
```

**Common misalignment patterns:**
- "X was significantly associated with Y" → but no P value provided
- "The model achieved high accuracy" → but AUC is only 0.72
- "We identified key genes" → but the selection criteria are unclear
- "Drug X reversed the disease signature" → but the enrichment score is marginal
- Figures referenced in text don't match (wrong panel letter)

### Lens 8: Logical Consistency

Check that no part of the manuscript contradicts another part.

**Cross-check:**
```
□ Abstract numbers match Results section numbers
□ Sample sizes in Methods match those in figure legends
□ Statistical methods described in Methods match those reported in Results
□ Conclusions in Discussion are supported by Results (no new claims)
□ Figure labels match text references (Fig. 1a in text = panel a in figure)
□ Gene/protein names are consistent throughout
□ p-value formatting is consistent
□ The number of features/genes/drugs is consistent across sections
```

### Lens 9: Impact & Significance Framing

Does the manuscript convincingly communicate WHY this work matters?

**Check each section for impact framing:**
- **Title**: Does it convey the advance? ("X reveals Y" not just "Analysis of X")
- **Abstract first sentence**: Does it hook a broad audience?
- **Introduction final paragraph**: Is the contribution clearly stated?
- **Discussion opening**: Is the main finding restated in broader context?
- **Discussion closing**: Does it look forward to future impact?

**The "so what?" test:**
After reading the whole paper, a reader should be able to answer:
1. What problem did they solve?
2. How did they solve it?
3. What did they find?
4. Why does it matter?
5. What should be done next?

If any of these is unclear, the framing needs work.

### Lens 10: Tone & Register

Is the writing appropriately formal and confident without being arrogant?

**Tone calibration:**
- **Too tentative**: "It might potentially be suggested that our results could possibly indicate..."
- **Too arrogant**: "We conclusively prove that..."
- **Just right**: "Our findings demonstrate that..." / "These results suggest that..."

**Register checks:**
- No colloquialisms ("a lot of", "kind of", "pretty significant")
- No hyperbole ("revolutionary", "groundbreaking", "unprecedented")
- No self-praise ("our novel and innovative approach")
- Appropriate confidence: match language strength to evidence strength
  - Strong evidence → "demonstrates", "establishes", "reveals"
  - Moderate evidence → "suggests", "indicates", "supports"
  - Weak evidence → "is consistent with", "may indicate"

### Lens 11: Completeness

Is anything missing that a reader would need?

**Missing content checklist:**
```
□ All abbreviations defined at first use
□ All statistical tests justified (why this test?)
□ All figure panels described in the text
□ All key prior work cited
□ Limitations honestly acknowledged
□ Alternative interpretations discussed
□ Data/code availability statements complete
□ Ethics statements included (if needed)
□ Funding acknowledged
□ All authors' contributions listed
```

### Lens 12: Reader Experience

The final lens — step back and evaluate the overall reading experience.

**Questions:**
- Is the paper enjoyable to read, or is it a chore?
- Does it maintain momentum throughout, or does it drag?
- Are there any sections where attention wanders?
- Is the length appropriate, or does it feel padded?
- Would you recommend this paper to a colleague?
- If you were a reviewer, would you accept this paper?

---

## Iteration Log Template

After each iteration, document findings:

```markdown
## Iteration [N]

**Date**: [YYYY-MM-DD]
**Focus**: [Which lenses revealed issues]

### Issues Found
1. [Section, paragraph] — [Description of issue] — [Severity: HIGH/MEDIUM/LOW]
2. ...

### Fixes Applied
1. [What was changed and why]
2. ...

### Remaining Concerns
- [Any issues that need further work]

### Assessment
- [ ] Ready for next iteration
- [ ] Ready for submission (no new concerns found)
```

---

## Convergence Criteria

The manuscript is ready for submission when ALL of the following are true:

```
□ A complete read-through through all 12 lenses produces ZERO new concerns
□ The naive reader test passes — a non-specialist can follow the entire paper
□ Every figure caption is self-explanatory without the main text
□ Every biological claim is supported by data AND literature
□ No AI vocabulary patterns remain
□ The abstract tells a complete story in ≤200 words
□ The 10-pass structural audit passes with no HIGH-severity issues
□ At least 3 complete iterations have been performed
□ The manuscript has been "read aloud" (mentally) for flow
□ You would be confident submitting this to Nature
```

---

## Common Patterns That Require Multiple Iterations

These issues are almost never caught in a single pass:

1. **Inconsistent terminology** — Using "disease signature", "gene set", "biomarker panel", and "feature set" interchangeably. Pick ONE term and use it throughout.

2. **Buried key findings** — The most important result is hidden in the middle of a paragraph instead of being the topic sentence.

3. **Missing "so what"** — Results are presented without interpretation. After every finding, the reader needs to know why it matters.

4. **Orphan references** — Citations that are mentioned once but never integrated into the argument.

5. **Figure-text mismatch** — The text describes a "clear separation" but the figure shows overlapping clusters.

6. **Scope creep in Discussion** — Starting to discuss topics not covered in the Results.

7. **Passive voice accumulation** — Individual passive sentences are fine, but three consecutive passive sentences kill momentum.

8. **Abstract drift** — The abstract no longer accurately reflects the paper after revisions to the main text.
