# Reviewer Expectations & Common Rejection Reasons

## How Nature Peer Review Works

### Editorial Triage (First 48 Hours)
- ~92% of submissions are desk-rejected at this stage
- Editors assess: scope, significance, accessibility, basic formatting
- The "90-second test": title + abstract + Figure 1 must convey the advance
- **Most common reason**: "Findings are incremental, not transformative"

### Peer Review (2–6 Weeks)
- Typically 2–3 reviewers selected by the editor
- Reviewers assess: scientific rigor, novelty, methodology, presentation
- Each reviewer provides: summary, major concerns, minor concerns, recommendation

### Decision Outcomes
| Decision | Meaning | Action |
|---|---|---|
| Accept | Rare on first submission | Celebrate |
| Minor revisions | Good science, small fixes needed | Fix and resubmit promptly |
| Major revisions | Promising but needs significant work | Address every concern thoroughly |
| Reject with resubmission | Not ready, but topic is interesting | Major overhaul, new data may be needed |
| Reject | Does not meet standards | Consider a different journal |

---

## Top 10 Reasons for Desk Rejection

Based on Nature editor feedback and editorial policies:

### 1. Insufficient Impact / Novelty
- **What they say**: "This is a solid contribution but may be more appropriate for a specialist journal."
- **Prevention**: Frame the significance broadly. Answer: "Why should a physicist care about this biology paper?"

### 2. Out of Scope
- **What they say**: "This work falls outside the scope of Nature."
- **Prevention**: Read 50 recent papers in Nature in your area. Does your work fit?

### 3. Inaccessible Writing
- **What they say**: "The manuscript is too specialized for Nature's broad readership."
- **Prevention**: First paragraph must be jargon-free. Use the funnel structure.

### 4. Incomplete Data
- **What they say**: "The findings would benefit from additional experiments/analyses."
- **Prevention**: Anticipate what validation experiments reviewers will request.

### 5. Methodological Concerns
- **What they say**: "The statistical/computational approach needs strengthening."
- **Prevention**: Use appropriate tests, report effect sizes, validate thoroughly.

### 6. Poor Presentation
- **What they say**: "The figures are difficult to interpret" or "The narrative is unclear."
- **Prevention**: Figures must be self-explanatory. Story arc must be clear.

### 7. Missing Controls
- **What they say**: "The study lacks appropriate controls or baselines."
- **Prevention**: Include negative controls, baseline comparisons, null models.

### 8. Overstated Claims
- **What they say**: "The conclusions are not fully supported by the data."
- **Prevention**: Match language strength to evidence strength. "Suggests" not "proves."

### 9. Ethical / Policy Issues
- **What they say**: "Missing data availability / ethics approval / competing interests."
- **Prevention**: Complete ALL required statements. Use the Pass 8 checklist.

### 10. Prior Publication / Preprint Issues
- **What they say**: "This work has been substantially published elsewhere."
- **Prevention**: Nature does accept papers with preprints, but disclose them.

---

## What Peer Reviewers Look For

### Reviewer Mindset
Reviewers are busy scientists who volunteer their time. They typically:
- Spend 2–4 hours reviewing a paper
- Read the abstract first, then figures, then methods
- Look for reasons to reject (it's faster than confirming validity)
- Focus on their area of expertise (they may miss things outside it)
- Write 1–2 pages of comments

### The Reviewer's Mental Checklist
```
1. Is this novel? Have I seen this before?
2. Is the methodology sound?
3. Are the statistics appropriate?
4. Are the conclusions supported by the data?
5. Are there alternative explanations?
6. Is it clearly written?
7. Are the figures clear and informative?
8. Would this paper change how I think or work?
```

### Reviewer Red Flags (What Triggers Negative Reviews)

| Red Flag | Reviewer Reaction |
|---|---|
| Missing error bars | "Authors don't understand variability" |
| No sample sizes in legends | "How can I assess statistical power?" |
| P < 0.05 without test name | "What test was used? Is it appropriate?" |
| "Data not shown" | "Why not? Are they hiding something?" |
| Claims without citations | "This is opinion, not evidence" |
| Figures with tiny labels | "Authors didn't check their own figures" |
| Obvious grammatical errors | "If they didn't proofread, what else did they rush?" |
| Discussion that repeats Results | "Author doesn't understand the purpose of Discussion" |
| No limitations section | "Author lacks self-awareness" |

---

## How to Write a Rebuttal (Response to Reviewers)

### Format
```
We thank the reviewers for their constructive comments, which have
significantly improved our manuscript. Below we address each point
in detail. Reviewer comments are in bold, our responses in regular text.

REVIEWER 1

**1. The authors should address the potential for overfitting given
the high dimensionality of the scRNA-seq data.**

We appreciate this important concern. To address overfitting, we have:
(1) Implemented five-fold stratified cross-validation...
(2) Used a held-out test set comprising 20% of cells...
(3) Added learning curves (new Extended Data Fig. X) showing...
The revised manuscript now includes these analyses (Results, page X,
lines Y–Z).

**2. The drug repurposing results need experimental validation.**

We agree that experimental validation would strengthen the findings.
However, as this is a computational study, in vitro/in vivo validation
is beyond its scope. We have:
(1) Acknowledged this as a limitation (Discussion, page X)...
(2) Proposed specific validation experiments for future work...
(3) Noted that 8 of our top candidates have published evidence...
```

### Rebuttal Best Practices
1. **Thank the reviewer** — even if the comment is harsh
2. **Address every point** — skipping one signals you can't answer it
3. **Be specific** — reference exact page/line numbers for changes
4. **Show, don't tell** — include new data, figures, or analyses
5. **Disagree respectfully** — "While we appreciate this perspective, our data suggest..."
6. **Never be defensive** — "As we clearly stated..." is adversarial
7. **Highlight major changes** — use a summary of changes at the top

### When Reviewers Are Wrong
Occasionally a reviewer misunderstands your work:
- **Clarify without condescension**: "We may not have been sufficiently clear. To address this..."
- **Provide evidence**: Cite literature supporting your approach
- **Improve the manuscript**: If a reviewer misunderstood, other readers might too
- **Never ignore**: Always respond, even if you disagree

---

## Pre-Submission Reviewer Simulation

Before submitting, role-play as three reviewers with different perspectives:

### Reviewer 1: The Methodologist
- "Are the methods appropriate and well-described?"
- "Is the statistical approach valid?"
- "Could I reproduce this work from the Methods section?"
- Focus areas: experimental design, sample sizes, controls, validation

### Reviewer 2: The Domain Expert
- "Is this novel in the context of prior literature?"
- "Are the biological interpretations correct?"
- "Are there important references missing?"
- Focus areas: literature context, biological plausibility, field impact

### Reviewer 3: The Skeptic
- "What's the weakest claim in this paper?"
- "What alternative explanation could account for these results?"
- "Is this result an artifact of the computational method?"
- Focus areas: confounders, artifacts, overstated conclusions

For each simulated reviewer, write 3 major concerns and address them in the manuscript BEFORE submission.
