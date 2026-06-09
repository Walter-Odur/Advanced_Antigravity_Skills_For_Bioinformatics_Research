"""
Manuscript Project Scaffolder
==============================
Creates a new manuscript project directory with all necessary files
and correct structure for the target journal.

Supports: LaTeX and DOCX output formats.

Usage:
    python scaffold_manuscript.py my_project --profile nature --format latex
    python scaffold_manuscript.py my_project --profile nature --format docx
    python scaffold_manuscript.py my_project --profile nature --format both
"""

import os
import sys
import argparse

os.chdir(os.path.dirname(os.path.abspath(__file__)))


LATEX_MAIN = r"""\documentclass[12pt]{article}

% ── Packages ──
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{mathptmx}          % Times New Roman
\usepackage[margin=1in]{geometry}
\usepackage{setspace}
\usepackage{graphicx}
\usepackage{xcolor}
\usepackage{hyperref}
\usepackage{natbib}
\usepackage{lineno}
\usepackage{amsmath}

% ── Settings ──
\doublespacing
\linenumbers
\hypersetup{colorlinks=true, linkcolor=blue, citecolor=blue, urlcolor=blue}

% ── Title ──
\title{[Your Title Here — ≤90 characters for Nature]}
\author{
    Author One\textsuperscript{1,*},
    Author Two\textsuperscript{1},
    Author Three\textsuperscript{2} \\[6pt]
    \small\textsuperscript{1}Department, Institution, City, Country \\
    \small\textsuperscript{2}Department, Institution, City, Country \\[6pt]
    \small\textsuperscript{*}Corresponding author: email@institution.edu
}
\date{}

\begin{document}
\maketitle

% ── Summary (Abstract) ──
% Nature: ≤200 words, single paragraph, bolded, NO references
\begin{abstract}
\textbf{%
[Context — one sentence establishing broad significance].
However, [specific knowledge gap].
Here we [approach/method], [using what data/tools].
We show that [main finding 1], and demonstrate that [main finding 2].
[Additional finding].
These results [implication], suggesting [broader impact].
}
\end{abstract}

% ── Main Body ──
% Nature Articles: ≤3,000 words (Introduction + Results + Discussion)
% Section headers are optional for Nature Articles

\section*{Introduction}
% ~500 words. Funnel: broad → narrow → gap → "Here we..."

[Hook: broad significance, accessible to all scientists.]

[Field context: prior work, 5–8 citations.]

[Gap: what remains unknown.]

Here we demonstrate that [your contribution].

\section*{Results}
% ~1,500 words. Each paragraph anchored to a figure.

\subsection*{[Descriptive subheading]}

[Result 1, referencing Fig. 1a.]

\subsection*{[Descriptive subheading]}

[Result 2, referencing Fig. 2.]

\section*{Discussion}
% ~1,000 words. Inverted funnel: narrow → broad.

[Restate main finding in broader context.]

[Compare to prior literature.]

[Mechanistic interpretation.]

[Limitations.]

[Future directions and implications.]

% ── Methods ──
% Nature: Online-only, no word limit, past tense
\section*{Methods}

\subsection*{Data sources}
[Describe data sources.]

\subsection*{Data processing}
[Describe preprocessing.]

\subsection*{Statistical analysis}
[Describe statistical methods, software versions.]

\subsection*{Code availability}
Custom code is available at [URL] and archived at Zenodo ([DOI]).

% ── References ──
\bibliographystyle{naturemag}
\bibliography{references}

% ── End Matter ──
\section*{Acknowledgements}
This work was supported by [funding source] (grant number [XXX]).

\section*{Author Contributions}
[A.O.] conceived the study. [B.C.] performed the analysis.
[A.O. and B.C.] wrote the manuscript. All authors reviewed the manuscript.

\section*{Competing Interests}
The authors declare no competing interests.

\section*{Data Availability}
[State where data can be accessed, with accession numbers/DOIs.]

% ── Figures ──
% Place figure environments here (legends at end for Nature)

\begin{figure}[htbp]
\centering
% \includegraphics[width=\textwidth]{figures/fig1.pdf}
\caption{\textbf{Fig. 1 $|$ [Title describing what the figure shows].}
\textbf{a}, [Panel a description]. \textbf{b}, [Panel b description].
$n$ = [sample size]. [Statistical test], [P value].
Error bars represent [s.d./s.e.m./95\% CI].}
\label{fig:1}
\end{figure}

\end{document}
"""

BIB_TEMPLATE = r"""% References for manuscript
% Use BibTeX format. See references/citation_formatting.md for Nature style.

@article{example2024,
  author  = {Last, First A and Other, Second B},
  title   = {Example article title},
  journal = {Nature},
  volume  = {600},
  pages   = {100--110},
  year    = {2024},
  doi     = {10.1038/example}
}
"""

DOCX_TEMPLATE_INFO = """
DOCX Manuscript Template
=========================
This directory is set up for a DOCX manuscript.

To create the manuscript:
1. Open Microsoft Word (or Google Docs)
2. Use the following style settings:
   - Font: 12-pt Times New Roman
   - Line spacing: Double
   - Margins: 1 inch (2.54 cm) all sides
   - Page numbers: Bottom center
   - Line numbers: Continuous (Layout → Line Numbers → Continuous)

3. Follow this section order:
   - Title (14-pt bold, centered)
   - Authors & Affiliations (12-pt, centered)
   - Summary (bold paragraph, ≤200 words)
   - Main Body (Introduction → Results → Discussion, ≤3,000 words)
   - Methods (no word limit)
   - References (numbered list, ≤50)
   - Acknowledgements
   - Author Contributions
   - Competing Interests
   - Data Availability
   - Code Availability
   - Figure Legends (grouped, each with bold title)

4. Save as .docx for submission

See references/nature_guide.md for complete formatting rules.
"""

GITIGNORE = """# LaTeX build files
*.aux
*.bbl
*.blg
*.log
*.out
*.toc
*.lof
*.lot
*.fls
*.fdb_latexmk
*.synctex.gz
*.dvi
*.ps

# OS files
.DS_Store
Thumbs.db

# IDE files
.vscode/
.idea/
"""


def scaffold(project_name, profile='nature', output_format='latex'):
    """Create a new manuscript project directory."""
    base = os.path.abspath(project_name)

    if os.path.exists(base):
        print(f"Error: Directory already exists: {base}")
        print("Choose a different name or delete the existing directory.")
        sys.exit(1)

    print(f"Creating manuscript project: {project_name}")
    print(f"Profile: {profile}")
    print(f"Format: {output_format}")
    print()

    # Create directories
    dirs = [
        base,
        os.path.join(base, 'figures'),
        os.path.join(base, 'supplementary'),
    ]
    if output_format in ('latex', 'both'):
        dirs.append(os.path.join(base, 'manuscript'))
    if output_format in ('docx', 'both'):
        dirs.append(os.path.join(base, 'manuscript_docx'))

    for d in dirs:
        os.makedirs(d, exist_ok=True)
        print(f"  📁 {os.path.relpath(d, os.path.dirname(base))}")

    # Create files
    if output_format in ('latex', 'both'):
        # main.tex
        tex_path = os.path.join(base, 'manuscript', 'main.tex')
        with open(tex_path, 'w', encoding='utf-8') as f:
            f.write(LATEX_MAIN)
        print(f"  📄 manuscript/main.tex")

        # references.bib
        bib_path = os.path.join(base, 'manuscript', 'references.bib')
        with open(bib_path, 'w', encoding='utf-8') as f:
            f.write(BIB_TEMPLATE)
        print(f"  📄 manuscript/references.bib")

    if output_format in ('docx', 'both'):
        # DOCX instructions
        docx_info = os.path.join(base, 'manuscript_docx', 'TEMPLATE_INFO.md')
        with open(docx_info, 'w', encoding='utf-8') as f:
            f.write(DOCX_TEMPLATE_INFO)
        print(f"  📄 manuscript_docx/TEMPLATE_INFO.md")

    # .gitignore
    gitignore_path = os.path.join(base, '.gitignore')
    with open(gitignore_path, 'w', encoding='utf-8') as f:
        f.write(GITIGNORE)
    print(f"  📄 .gitignore")

    # README
    readme_path = os.path.join(base, 'README.md')
    with open(readme_path, 'w', encoding='utf-8') as f:
        f.write(f"# {project_name}\n\n")
        f.write(f"Scientific manuscript project targeting **{profile.title()}**.\n\n")
        f.write("## Structure\n\n")
        f.write("```\n")
        f.write(f"{project_name}/\n")
        if output_format in ('latex', 'both'):
            f.write("├── manuscript/\n")
            f.write("│   ├── main.tex          # Main LaTeX document\n")
            f.write("│   └── references.bib    # BibTeX references\n")
        if output_format in ('docx', 'both'):
            f.write("├── manuscript_docx/\n")
            f.write("│   └── TEMPLATE_INFO.md  # DOCX formatting guide\n")
        f.write("├── figures/              # High-resolution figures\n")
        f.write("├── supplementary/        # Supplementary materials\n")
        f.write("└── .gitignore\n")
        f.write("```\n\n")
        f.write("## Quality Check\n\n")
        f.write("```bash\n")
        if output_format in ('latex', 'both'):
            f.write(f"python quality_check.py manuscript/main.tex --profile {profile}\n")
            f.write(f"python word_counter.py manuscript/main.tex --profile {profile}\n")
            f.write(f"python citation_validator.py manuscript/references.bib manuscript/main.tex\n")
        f.write("```\n")
    print(f"  📄 README.md")

    print(f"\n✅ Project scaffolded successfully: {base}")
    print(f"\nNext steps:")
    if output_format in ('latex', 'both'):
        print(f"  1. Edit manuscript/main.tex")
        print(f"  2. Add references to manuscript/references.bib")
    if output_format in ('docx', 'both'):
        print(f"  1. Create DOCX following manuscript_docx/TEMPLATE_INFO.md")
    print(f"  3. Add figures to figures/")
    print(f"  4. Run quality_check.py when ready")


def main():
    parser = argparse.ArgumentParser(
        description='Scaffold a new manuscript project'
    )
    parser.add_argument('project', help='Project directory name')
    parser.add_argument('--profile', '-p', default='nature',
                        help='Journal profile (default: nature)')
    parser.add_argument('--format', '-f', default='both',
                        choices=['latex', 'docx', 'both'],
                        help='Output format (default: both)')

    args = parser.parse_args()
    scaffold(args.project, args.profile, args.format)


if __name__ == '__main__':
    main()
