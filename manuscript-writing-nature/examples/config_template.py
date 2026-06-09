"""
Manuscript Writing — Configuration Template
=============================================
Copy this file to your manuscript project and customize.
All parameters here correspond to values used in the scripts.
"""

# ── Journal Profile ──
JOURNAL = "nature"                    # Options: nature (more coming)

# ── Output Format ──
OUTPUT_FORMAT = "both"                # Options: latex, docx, both

# ── Word Limits (Nature defaults) ──
TITLE_MAX_CHARS = 90
SUMMARY_MAX_WORDS = 200
MAIN_BODY_MAX_WORDS = 3000            # Intro + Results + Discussion
LEGEND_MAX_WORDS = 300                # Per figure legend
MAX_REFERENCES = 50
MAX_DISPLAY_ITEMS = 6                 # Figures + tables combined
MAX_EXTENDED_DATA = 10

# ── File Paths ──
MANUSCRIPT_DIR = "manuscript"
FIGURES_DIR = "figures"
SUPPLEMENTARY_DIR = "supplementary"
MAIN_TEX = "manuscript/main.tex"
REFERENCES_BIB = "manuscript/references.bib"

# ── Figure Settings ──
FIGURE_DPI = 600                      # For publication quality
SINGLE_COLUMN_MM = 89                 # Nature single column
DOUBLE_COLUMN_MM = 183                # Nature double column
FULL_PAGE_HEIGHT_MM = 247
FIGURE_FONT = "Arial"
FIGURE_FONT_SIZE = 8                  # Points
PANEL_LABEL_SIZE = 10                 # Points

# ── Color Palette (Color-blind safe, Wong 2011) ──
COLORS = {
    'blue': '#0072B2',
    'orange': '#E69F00',
    'green': '#009E73',
    'vermillion': '#D55E00',
    'sky_blue': '#56B4E9',
    'yellow': '#F0E442',
    'black': '#000000',
    'grey': '#999999',
}

# ── Quality Check Settings ──
MAX_SENTENCE_WORDS = 40               # Flag sentences longer than this
MAX_PASSIVE_VOICE_PCT = 30            # Warn if passive voice exceeds this %
SELF_CITATION_THRESHOLD = 20          # Warn if self-citations exceed this %

# ── LaTeX Compilation ──
LATEX_ENGINE = "pdflatex"             # Options: pdflatex, xelatex, lualatex
BIBTEX_ENGINE = "bibtex"              # Options: bibtex, biber
