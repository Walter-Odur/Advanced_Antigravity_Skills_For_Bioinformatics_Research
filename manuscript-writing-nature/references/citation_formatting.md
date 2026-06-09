# Citation Formatting for Nature

## In-Text Citation Style

Nature uses **numbered superscript** citations:

```latex
... as demonstrated previously\textsuperscript{1}.
... multiple studies\textsuperscript{1--3} have confirmed...
... prior work\textsuperscript{1,4,7} suggests...
```

### Rules
- Numbers are **superscript**, placed **after punctuation** (period, comma)
- Consecutive references: use en-dash (1–3), not (1,2,3)
- Non-consecutive: separate with commas (1,4,7)
- Numbers assigned **in order of first appearance** in the text
- Maximum 50 references for Nature Articles

---

## Reference List Format

### Journal Article (≤5 authors — list all)
```
1. Subramanian, A., Narayan, R., Corsello, S. M., Peck, D. D.
   & Bhatt, D. L. A next generation connectivity map: L1000
   platform and the first 1,000,000 profiles. Cell 171,
   1437–1452 (2017).
```

### Journal Article (6+ authors — first author et al.)
```
2. Heumos, L. et al. Best practices for single-cell analysis
   across modalities. Nat. Rev. Genet. 24, 550–572 (2023).
```

### Book
```
3. Hastie, T., Tibshirani, R. & Friedman, J. The Elements of
   Statistical Learning 2nd edn (Springer, 2009).
```

### Book Chapter
```
4. Smith, J. in Computational Biology (ed. Jones, K.) 45–78
   (Academic Press, 2020).
```

### Preprint
```
5. Doe, J. et al. Novel approach to drug repurposing. Preprint
   at https://doi.org/10.1101/2024.01.15.12345 (2024).
```

### Conference Paper
```
6. Chen, T. & Guestrin, C. XGBoost: a scalable tree boosting
   system. In Proc. 22nd ACM SIGKDD International Conference on
   Knowledge Discovery and Data Mining 785–794 (ACM, 2016).
```

### Software
```
7. Pedregosa, F. et al. Scikit-learn: machine learning in Python.
   J. Mach. Learn. Res. 12, 2825–2830 (2011).
```

---

## Formatting Rules

### Author Names
- **Format**: Last, First Initial. Middle Initial.
- **Separator**: Commas between authors, **&** before last author
- **6+ authors**: First author followed by "et al."
- **Consistency**: Use the same name format throughout

### Journal Names
- **Italicized** and **abbreviated** using MEDLINE standard
- Common abbreviations:

| Full Name | Abbreviation |
|---|---|
| Nature | *Nature* |
| Science | *Science* |
| Cell | *Cell* |
| Nature Methods | *Nat. Methods* |
| Nature Machine Intelligence | *Nat. Mach. Intell.* |
| Nature Reviews Genetics | *Nat. Rev. Genet.* |
| Nature Reviews Drug Discovery | *Nat. Rev. Drug Discov.* |
| Nature Communications | *Nat. Commun.* |
| Nature Biotechnology | *Nat. Biotechnol.* |
| Genome Biology | *Genome Biol.* |
| Genome Research | *Genome Res.* |
| Molecular Systems Biology | *Mol. Syst. Biol.* |
| Bioinformatics | *Bioinformatics* |
| BMC Genomics | *BMC Genomics* |
| PNAS | *Proc. Natl Acad. Sci. USA* |
| JAMA | *JAMA* |
| The Lancet | *Lancet* |
| New England Journal of Medicine | *N. Engl. J. Med.* |
| Journal of Machine Learning Research | *J. Mach. Learn. Res.* |
| Physical Review Letters | *Phys. Rev. Lett.* |
| Annual Review of Immunology | *Annu. Rev. Immunol.* |
| Machine Learning | *Mach. Learn.* |
| Nucleic Acids Research | *Nucleic Acids Res.* |

### Volume, Pages, Year
- **Volume**: In bold (in print), regular in references: just the number
- **Pages**: Start–end with en-dash (not hyphen): 1437–1452
- **Year**: In parentheses at the end: (2017).
- **Article number**: For journals without page numbers: use article number
  - Example: `Nat. Commun. 12, 4567 (2021).`

### DOIs
- Nature references typically do NOT include DOIs in the reference list
- DOIs are used internally for linking but not printed
- Exception: preprints MUST include DOI

---

## BibTeX for Nature

### BibTeX Style File
Use `naturemag.bst` or a compatible style:

```latex
\bibliographystyle{naturemag}
\bibliography{references}
```

### BibTeX Entry Examples

```bibtex
@article{subramanian2017,
  author  = {Subramanian, Aravind and Narayan, Rajiv and Corsello,
             Steven M and Peck, David D and Natoli, Ted E and Lu,
             Xiaodong and Gould, Joshua and Davis, John F and
             Tubelli, Andrew A and Asiedu, Jacob K and others},
  title   = {A next generation connectivity map: {L1000} platform
             and the first 1,000,000 profiles},
  journal = {Cell},
  volume  = {171},
  number  = {6},
  pages   = {1437--1452},
  year    = {2017},
  doi     = {10.1016/j.cell.2017.10.049}
}

@article{heumos2023,
  author  = {Heumos, Lukas and Schaar, Anna C and Lance,
             Christopher and others},
  title   = {Best practices for single-cell analysis across
             modalities},
  journal = {Nature Reviews Genetics},
  volume  = {24},
  pages   = {550--572},
  year    = {2023},
  doi     = {10.1038/s41576-023-00586-w}
}

@book{hastie2009,
  author    = {Hastie, Trevor and Tibshirani, Robert and
               Friedman, Jerome},
  title     = {The Elements of Statistical Learning},
  edition   = {2nd},
  publisher = {Springer},
  year      = {2009}
}

@inproceedings{chen2016xgboost,
  author    = {Chen, Tianqi and Guestrin, Carlos},
  title     = {{XGBoost}: A Scalable Tree Boosting System},
  booktitle = {Proceedings of the 22nd ACM SIGKDD International
               Conference on Knowledge Discovery and Data Mining},
  pages     = {785--794},
  year      = {2016},
  publisher = {ACM}
}
```

### Protecting Capitalization
Use braces `{}` around words that must remain capitalized:
```bibtex
title = {{XGBoost}: A Scalable Tree Boosting System}
title = {A next generation connectivity map: {L1000} platform}
title = {{SHAP} values for interpretable machine learning}
```

---

## Common Citation Mistakes

| Mistake | Correction |
|---|---|
| Citing reviews for specific findings | Cite the original research paper |
| Only citing recent work | Include foundational/seminal papers |
| Excessive self-citation (>20%) | Ensure balanced citation profile |
| Missing citations for methods | Cite software/tools in Methods |
| Citing retracted papers | Check Retraction Watch database |
| Wrong year/volume/pages | Verify against DOI lookup |
| Inconsistent format | Use BibTeX with consistent .bst file |

## Citation Density Guidelines

| Section | Typical Citations |
|---|---|
| Introduction | 15–20 references (context, prior work, gap) |
| Results | 3–5 references (comparisons, benchmarks) |
| Discussion | 15–20 references (comparison to literature) |
| Methods | 10–15 references (tools, software, protocols) |
| **Total** | **≤50 for Nature Articles** |
