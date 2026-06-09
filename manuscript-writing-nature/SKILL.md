---
name: manuscript-writing
description: "Advanced scientific manuscript writing skill with journal-specific profiles, narrative writing guidance, an automated 10-pass structural quality audit, and a deep iterative AI self-review engine with 12 critical lenses. Supports LaTeX and DOCX output. Currently supports Nature; designed to expand to Cell, Science, Lancet, NEJM, and more."
metadata:
    skill-author: Walter Odur
    version: "1.0.0"
    domain: scientific-writing, academic-publishing, manuscript-preparation
risk: low
source: custom
---

# Scientific Manuscript Writing Skill

## Overview

A comprehensive, publication-grade manuscript writing system that guides the creation of scientific papers from outline to submission-ready manuscripts. Each manuscript is built using **journal-specific profiles** that encode exact formatting rules, word limits, writing style requirements, and submission checklists from official author guidelines.

The skill implements a **two-tier quality assurance system**:
1. **10-pass structural quality audit** — automated and manual checks for formatting, word counts, citations, and compliance
2. **Deep iterative AI self-review** — unlimited-iteration critical reading through 12 lenses (narrative flow, figure captions, biological interpretation, AI vocabulary detection, and more) until no new concerns remain

### Currently Supported Journals

| Journal | Profile | Status |
|---|---|---|
| **Nature** | `journal_profiles/nature.yaml` | ✅ Complete |
| Cell | — | 🔜 Planned |
| Science | — | 🔜 Planned |
| Lancet | — | 🔜 Planned |
| NEJM | — | 🔜 Planned |
| Nature Methods | — | 🔜 Planned |
| Genome Biology | — | 🔜 Planned |

## When to Use This Skill

Activate this skill when the user needs to:

- **Write** a new manuscript for a specific journal (Nature, Cell, etc.)
- **Revise** an existing manuscript to meet journal requirements
- **Quality-check** a manuscript for compliance and readability
- **Format** a manuscript (LaTeX/BibTeX or DOCX) for submission
- **Scaffold** a new manuscript project with correct structure (`scripts/scaffold_manuscript.py`)
- **Prepare** submission materials (cover letter, checklists, statements)
- **Respond** to reviewer comments with a structured rebuttal
- **Convert** a manuscript between journal formats

## Core Workflow: The 8-Step Manuscript Pipeline

> **IMPORTANT**: Before each step, read the relevant reference document listed in that step's instructions. Do not skip reference documents — they contain critical journal-specific requirements.

### Step 1: Journal Selection & Scope Check

Before writing a single word, verify journal fit:

```
CHECKLIST:
□ Does the research represent a significant advance?
□ Will it interest readers beyond the immediate subfield?
□ Does the methodology meet the journal's rigor standards?
□ Is the scope appropriate (basic science, clinical, computational)?
□ Have similar papers been published in this journal recently?
```

**For Nature specifically:**
- Acceptance rate is <8% — only truly transformative findings
- Must have broad interdisciplinary significance
- Editors triage in <48 hours based on title + abstract + Figure 1
- The "90-second test": Can a non-specialist editor grasp the significance in 90 seconds?

### Step 2: Build the Story Arc

Every manuscript needs a compelling narrative. Use this framework:

```
THE STORY ARC:
1. HOOK        → First sentence captures the reader ("Why should I care?")
2. CONTEXT     → Broader field significance (2-3 sentences)
3. GAP         → What we don't know / what's missing (1-2 sentences)
4. APPROACH    → What we did to address the gap (1-2 sentences)
5. KEY FINDING → The headline result (1-2 sentences)
6. IMPACT      → Why this matters for the field and beyond (1-2 sentences)
```

This arc applies at EVERY level:
- **Whole paper**: Abstract follows this exact arc
- **Each section**: Introduction builds hook→gap, Results delivers findings
- **Each paragraph**: Topic sentence → evidence → interpretation
- **Each figure**: Caption tells a mini-story

### Step 3: Figures First

**Build your visual narrative before writing text.** This is the most critical but least-followed principle:

1. Draft all figures with complete captions
2. Arrange figures in logical story order
3. Write figure legends that are self-explanatory
4. The figures alone (without text) should tell the complete story
5. Only THEN write the text to connect and explain the figures

**Nature Figure Rules (from `journal_profiles/nature.yaml`):**
- Maximum 6 display items (figures + tables combined)
- Extended Data: up to 10 additional items
- Single column: 89mm; Double column: 183mm
- Minimum 300 DPI (600 for line art)
- Arial/Helvetica labels, lowercase panel letters (a, b, c)
- Color-blind safe palettes required

See `references/figure_standards.md` for complete specifications.

### Step 4: Section-by-Section Drafting

**Write sections in this order** (NOT the order they appear in the paper):

```
DRAFTING ORDER:
1. Methods     → Easiest to write; establishes what you did
2. Results     → Describe findings objectively, reference figures
3. Discussion  → Interpret results, compare to literature, limitations
4. Introduction→ Now you know what story to set up
5. Abstract    → Compress the entire paper into 200 words
6. Title       → Distill the key message into <90 characters
```

**For each section, follow the two-stage process:**

**Stage 1: Outline** — Create bullet points with key claims, data, and citations
**Stage 2: Prose** — Convert every bullet point into flowing paragraphs

> **CRITICAL: Never leave bullet points in the final manuscript.** All sections must be written as complete, flowing paragraphs with proper transitions.

#### Section Guidelines

**Title** (see `references/nature_guide.md`):
- ≤90 characters for Nature
- Active and specific: "X reveals Y" or "X enables Y"
- No abbreviations, formulae, or trade names
- Must stand alone without the abstract

**Summary/Abstract** (see `references/nature_guide.md`):
- Single unstructured paragraph, ≤200 words for Nature
- **Bolded** in Nature format
- NO references, NO labels (Background:, Methods:, etc.)
- Follows the story arc: context → gap → approach → findings → impact
- Every word must earn its place — this is the most-read section

**Introduction** (see `references/writing_style.md`):
- ~500 words for Nature Articles
- Paragraph 1: Hook + broad context (accessible to ALL scientists)
- Paragraph 2: Narrow to specific field, cite key prior work
- Paragraph 3: Identify the gap, state what remains unknown
- Final sentence: "Here we show/report/demonstrate..." (your approach)
- NO results here — save them for the Results section

**Results** (see `references/nature_guide.md`):
- ~1,500 words for Nature Articles
- Each paragraph: one logical unit tied to a figure panel
- Pattern: "We found that X (Fig. 1a). This result suggests Y because Z."
- Past tense for your findings
- Brief interpretation is acceptable ("suggesting that...") but save deep interpretation for Discussion
- Reference every figure panel explicitly

**Discussion** (see `references/writing_style.md`):
- ~1,000 words for Nature Articles
- Paragraph 1: Restate the main finding in broader context
- Paragraph 2-3: Compare to prior literature, explain discrepancies
- Paragraph 4: Mechanistic interpretation or model
- Paragraph 5: Limitations (acknowledge honestly — reviewers will find them)
- Final paragraph: Future directions and broader implications
- Do NOT introduce new data
- Do NOT repeat Results verbatim

**Methods** (see `references/nature_guide.md`):
- Online-only for Nature (no word limit)
- Written in past tense
- Enough detail for reproducibility
- Statistical methods: test names, software, thresholds, corrections
- Include subsections: Study design, Data processing, ML models, Statistical analysis
- Ethics statements if applicable

**Required Statements** (see `references/data_code_availability.md`):
- Data Availability Statement (mandatory)
- Code Availability Statement (mandatory)
- Author Contributions (CRediT taxonomy or free text)
- Competing Interests Declaration
- Ethics Approval (if applicable)

### Step 5a: The 10-Pass Structural Quality Audit

After drafting, run 10 structural quality passes. Each pass has a specific focus. Fix ALL issues found before moving to the next pass.

See `references/quality_checklist.md` for the complete protocol.

**Quick Reference:**

| Pass | Focus | Tool |
|---|---|---|
| 1 | Structure & completeness | Manual |
| 2 | Word count compliance | `scripts/word_counter.py` |
| 3 | Narrative flow & story arc | Manual |
| 4 | Accessibility & jargon | Manual |
| 5 | Figure/table compliance | Manual |
| 6 | Statistical reporting | Manual + `references/statistical_reporting.md` |
| 7 | Citation integrity | `scripts/citation_validator.py` |
| 8 | Required statements | Manual + checklist |
| 9 | Consistency (terms, tense, format) | Manual |
| 10 | Reviewer simulation | Manual |

**Run the automated quality check (LaTeX or DOCX):**
```bash
python scripts/quality_check.py manuscript/main.tex --profile nature
python scripts/quality_check.py manuscript/paper.docx --profile nature
```

**Scaffold a new project (LaTeX, DOCX, or both):**
```bash
python scripts/scaffold_manuscript.py my_project --profile nature --format both
```

### Step 5b: Deep Iterative AI Self-Review (MANDATORY)

**This is the skill's most critical innovation.** After the structural audit passes, perform a deep iterative self-review that goes far beyond formatting checks. This review examines the manuscript through 12 critical lenses and repeats until NO new concerns are found.

See `references/ai_review_protocol.md` for the complete protocol.

**The fundamental question on each iteration:**

> *"Critically and deeply read through the final manuscript. Assuming you are an absolutely new reader with no idea how the analysis was performed — would you understand easily, without any difficulty?"*

**The 12 Review Lenses (applied on EVERY iteration):**

| # | Lens | What to Check |
|---|---|---|
| 1 | Narrative flow & coherence | Logical transitions, story momentum, no "orphan" paragraphs |
| 2 | Figure captions | Self-explanatory without main text? n, tests, p-values included? |
| 3 | Biological interpretation | Grounded in high-impact literature? Mechanistically plausible? |
| 4 | AI vocabulary detection | Eliminate "delve into", "leveraging", "tapestry", "multifaceted" |
| 5 | Clarity & precision | Unambiguous pronouns, specific numbers, no run-on sentences |
| 6 | Methods reproducibility | Could someone reproduce the analysis from Methods alone? |
| 7 | Data-claim alignment | Every claim traceable to a figure, statistic, or citation |
| 8 | Logical consistency | Abstract numbers match Results, Methods match figure legends |
| 9 | Impact framing | "So what?" test — significance clear in every section |
| 10 | Tone & register | Confident but not arrogant, no hyperbole, appropriate hedging |
| 11 | Completeness | All abbreviations defined, limitations acknowledged, alternatives discussed |
| 12 | Reader experience | Enjoyable to read? Maintains momentum? Would you accept as reviewer? |

**Iteration protocol:**
```
REPEAT {
    1. Read the ENTIRE manuscript from title to references
    2. Apply ALL 12 review lenses simultaneously
    3. Document every concern found
    4. Fix ALL concerns
    5. Re-read fixed sections in context
} UNTIL a full read-through produces ZERO new concerns
```

**Minimum**: 3 complete iterations (even if clean)
**Exit condition**: Full read-through with ZERO new concerns

**Study real Nature papers** as the quality benchmark — see `references/exemplary_articles.md` for curated examples from AlphaFold, CMap L1000, Scanpy, and other landmark papers.

### Step 6: Formatting & Compliance

After the quality audit (Steps 5a + 5b):

1. **LaTeX compilation** — verify clean build with no warnings
2. **Reference formatting** — check against Nature style
3. **Figure quality** — verify DPI, dimensions, color-blind safety
4. **Line numbers** — required for Nature submission
5. **Double spacing** — required for initial submission
6. **Page limits** — verify total length

**Supplementary Information:**
- Move large datasets, additional validations, and secondary analyses to Supplementary
- Reference as "Supplementary Fig. 1", "Supplementary Table 1"
- Supplementary Information is NOT peer-reviewed to the same standard as Extended Data
- Extended Data (≤10 items) IS peer-reviewed — use for important supporting results

### Step 7: Cover Letter & Submission

Draft a cover letter that:
- Addresses the editor by name (if known)
- States the article type and title
- Summarizes the significance in 2-3 sentences (non-technical)
- Explains why this paper fits the journal's scope
- Lists potential reviewers (and any to exclude)
- Confirms compliance with all policies

**Submission Checklist:**
```
□ Manuscript file (PDF or LaTeX)
□ Cover letter
□ All figure files (high-resolution)
□ Supplementary Information (if any)
□ Nature Reporting Summary (completed form)
□ Data Availability Statement in manuscript
□ Code Availability Statement in manuscript
□ Author contributions in manuscript
□ Competing interests declaration
□ Ethics approval documentation (if applicable)
□ ORCID IDs for all authors
```

### Step 8: Handling Revisions

If the manuscript receives a revision decision (major or minor):

1. **Read ALL reviewer comments** before responding to any single one
2. **Create a point-by-point response document** (see `references/reviewer_expectations.md`)
3. **Address every single comment** — even if you disagree, explain why
4. **Track changes** in the manuscript (LaTeX: use `latexdiff`; DOCX: use Track Changes)
5. **Provide both tracked and clean versions**
6. **Update figures** if any data presentation was criticized
7. **Re-run the full quality audit** (Steps 5a + 5b) on the revised manuscript
8. **Update Data/Code Availability** if new analyses were added

> **CRITICAL**: The rebuttal letter is as important as the manuscript itself. See `references/reviewer_expectations.md` for format, tone, and strategy.

---

## The Quality-First Philosophy

### Why 10 Passes?

Most rejected manuscripts fail not because of bad science, but because of **presentation issues** that are detectable before submission:

| Rejection Reason | Detection Pass |
|---|---|
| Findings not significant enough | Pass 3 (narrative), Pass 10 (reviewer sim) |
| Poor writing quality | Pass 4 (accessibility), Pass 9 (consistency) |
| Missing/incomplete statistics | Pass 6 (statistics) |
| Missing required statements | Pass 8 (statements) |
| Figures unclear or non-compliant | Pass 5 (figures) |
| Over word limit | Pass 2 (word count) |
| Insufficient/incorrect citations | Pass 7 (citations) |

### The "Read-Through" Test

After completing all quality passes (Steps 5a + 5b), perform the final "read-through" test:

1. **Read the entire manuscript aloud** — Does it flow naturally?
2. **Read ONLY the abstract** — Does it tell a complete, compelling story?
3. **Read ONLY the figure legends** — Do the figures alone tell the story?
4. **Read ONLY the first sentence of each paragraph** — Does the logic flow?
5. **Give it to a colleague outside your field** — Can they understand the significance?

---

## Integration with Other Skills

This skill works with:

| Skill | Integration |
|---|---|
| `scientific-writing` | General IMRAD principles, writing style fundamentals |
| `citation-management` | BibTeX creation, validation, DOI-to-reference conversion |
| `scrna-drug-repurposing` | Source data/results for the TB manuscript example |
| `matplotlib` / `seaborn` | Figure generation meeting journal specs |
| `latex-paper-conversion` | Converting between journal formats |

---

## Reference Documentation

| Document | Purpose |
|---|---|
| `references/nature_guide.md` | Complete Nature formatting rules and requirements |
| `references/writing_style.md` | Nature narrative writing: hooks, accessibility, storytelling |
| `references/quality_checklist.md` | Detailed 10-pass structural audit protocol |
| `references/ai_review_protocol.md` | **Deep iterative AI self-review with 12 lenses** |
| `references/exemplary_articles.md` | **Real top Nature/Cell/Science papers as writing models** |
| `references/statistical_reporting.md` | Statistical reporting standards for Nature |
| `references/data_code_availability.md` | FAIR principles, repository requirements |
| `references/reviewer_expectations.md` | Common rejection reasons, reviewer psychology |
| `references/figure_standards.md` | Figure specifications across journals |
| `references/citation_formatting.md` | Nature citation style and BibTeX formatting |

## Scripts

| Script | Purpose |
|---|---|
| `scripts/quality_check.py` | Automated 10-pass quality audit engine |
| `scripts/scaffold_manuscript.py` | Scaffold a new manuscript project |
| `scripts/word_counter.py` | Section-by-section word counting for LaTeX |
| `scripts/citation_validator.py` | Validate citations against .bib file |

## Journal Profiles

| Profile | Status |
|---|---|
| `journal_profiles/nature.yaml` | ✅ Complete |

## Examples

| Example | Description |
|---|---|
| `examples/tb_manuscript/` | Template structure for a Nature Article |
| `examples/config_template.py` | Configuration template for new manuscripts |

---

## Limitations

- This skill provides writing guidance and quality checking — it does not replace peer review
- Journal guidelines change periodically; verify current requirements before submission
- The automated quality checker catches formatting issues; the AI self-review addresses narrative quality
- Biological interpretation must be validated by domain experts
- The deep iterative review converges faster with user feedback between iterations
