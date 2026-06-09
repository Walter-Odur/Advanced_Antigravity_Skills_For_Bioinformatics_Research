# Figure Standards for Nature

## Dimensions

| Layout | Width | Use Case |
|---|---|---|
| Single column | 89 mm (3.50 in / 252 pt) | Simple panels, single graphs |
| 1.5 columns | 120 mm (4.72 in / 340 pt) | Medium complexity |
| Double column | 183 mm (7.20 in / 519 pt) | Multi-panel figures |
| Full page height | 247 mm (9.72 in / 700 pt) | Maximum figure height |

## Resolution

| Content Type | Minimum DPI | Recommended DPI |
|---|---|---|
| Photographs, micrographs | 300 | 600 |
| Line art, vector graphics | 600 | 1200 |
| Combination (line + photo) | 500 | 600 |

## File Formats

| Stage | Accepted Formats | Preferred |
|---|---|---|
| Initial submission | JPEG, PNG, TIFF, PDF, EPS | PDF or TIFF |
| Revision / Accepted | TIFF, EPS, PDF | TIFF (lossless) |

**Never use**: BMP, GIF, PowerPoint slides as figures

## Typography in Figures

- **Font**: Arial or Helvetica (sans-serif only)
- **Minimum font size**: 5 pt (after scaling to final size)
- **Recommended size**: 7–10 pt for labels, 8–12 pt for panel labels
- **Panel labels**: Lowercase bold (**a**, **b**, **c**, ...), top-left of each panel
- **Axis labels**: Include units in parentheses: "Expression level (log₂ CPM)"
- **Consistent sizing**: Same font size across ALL panels in ALL figures

## Color Standards

### Color-Blind Safe Palettes

Approximately 8% of males have some form of color vision deficiency. Use palettes distinguishable by all viewers.

**Recommended palette** (Wong, Nature Methods 2011):
| Color | Hex | RGB | Use |
|---|---|---|---|
| Blue | #0072B2 | (0, 114, 178) | Primary |
| Orange | #E69F00 | (230, 159, 0) | Secondary |
| Green | #009E73 | (0, 158, 115) | Tertiary |
| Vermillion | #D55E00 | (213, 94, 0) | Alert/highlight |
| Sky blue | #56B4E9 | (86, 180, 233) | Background/light |
| Yellow | #F0E442 | (240, 228, 66) | Accent (use sparingly) |
| Black | #000000 | (0, 0, 0) | Text/borders |
| Grey | #999999 | (153, 153, 153) | Neutral/background |

### Rules
- **Never** use red-green as the only distinguishing feature
- **Always** add shape/pattern differences alongside color
- Use **sequential** colormaps (viridis, plasma) for continuous data
- Use **diverging** colormaps (RdBu_r) for data centered on zero
- Test with a color-blindness simulator (e.g., Color Oracle)

## Multi-Panel Figure Layout

### Panel Organization
```
┌─────────┬─────────┐
│    a    │    b    │   Top row: overview / primary results
├─────────┼─────────┤
│    c    │    d    │   Bottom row: supporting / detailed results
└─────────┴─────────┘
```

### Best Practices
- Logical flow: left→right, top→bottom
- Related panels adjacent to each other
- Consistent axis scales across comparable panels
- Shared legends when possible (reduces clutter)
- White space between panels (not too crowded)
- Panel labels (**a**, **b**) consistently placed

## Figure Types & When to Use Them

| Data Type | Recommended Plot | When to Use |
|---|---|---|
| Distribution (1 group) | Histogram, violin | Showing spread of values |
| Distribution (2+ groups) | Box plot, violin, bee swarm | Comparing distributions |
| Comparison (categories) | Bar chart (with data points) | Comparing means across groups |
| Trend over time | Line graph | Time series, dose-response |
| Correlation (2 variables) | Scatter plot | Relationships between variables |
| High-dimensional | Heatmap | Gene expression, correlation matrices |
| Dimensionality reduction | UMAP, t-SNE, PCA scatter | Clustering, population structure |
| Model performance | ROC curves, PR curves | Classification evaluation |
| Feature importance | Horizontal bar chart | SHAP values, ranked features |
| Composition | Stacked bar, donut chart | Proportions, cell type composition |
| Network | Network graph | Protein interactions, pathways |

### Modern Best Practices
- **Show individual data points** on bar charts (not just the bar)
- **Use violin plots** instead of bar charts for distributions
- **Avoid pie charts** — use bar charts or stacked bars instead
- **Avoid 3D plots** unless the third dimension adds information
- **Use consistent colors** across all figures for the same groups
- **Avoid chart junk**: gridlines, shadows, 3D effects, gradients

## Figure Legend Structure

Every legend must follow this structure:

```
**Fig. X | [One-sentence title describing the panel's content].**
**a**, [Description of panel a, including what is plotted, what
colors/shapes represent, and key observation]. **b**, [Description
of panel b]. n = [sample size] [units]. Statistical test: [name],
[correction if applicable]. Error bars represent [s.d. / s.e.m. /
95% CI]. [Scale bar information if applicable].
```

### Legend Word Limits
| Journal | Limit per Legend |
|---|---|
| Nature | ≤300 words |
| Nature Methods | ≤300 words |
| Cell | ≤250 words |
| Science | ≤200 words |

## Generating Publication-Quality Figures

### Python (Matplotlib/Seaborn)

```python
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Nature-compliant settings
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 8,
    'axes.labelsize': 9,
    'axes.titlesize': 10,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'legend.fontsize': 7,
    'figure.dpi': 300,
    'savefig.dpi': 600,
    'savefig.bbox': 'tight',
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.linewidth': 0.5,
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
})

# Single column figure (89mm = 3.5 inches)
fig, ax = plt.subplots(figsize=(3.5, 2.5))

# Double column figure (183mm = 7.2 inches)
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))

# Save in multiple formats
plt.savefig('fig1.pdf', dpi=600, facecolor='white')
plt.savefig('fig1.tiff', dpi=600, facecolor='white', pil_kwargs={'compression': 'tiff_lzw'})
plt.savefig('fig1.png', dpi=600, facecolor='white')
```

### R (ggplot2)
```r
library(ggplot2)

# Nature-compliant theme
theme_nature <- function() {
  theme_minimal() +
  theme(
    text = element_text(family = "Arial", size = 8),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 7),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2, color = "#EEEEEE")
  )
}

# Save at correct dimensions
ggsave("fig1.pdf", width = 89, height = 60, units = "mm", dpi = 600)
```

## Extended Data vs. Supplementary Information

Nature distinguishes between two types of supporting display items:

### Extended Data (Peer-Reviewed)
- **Maximum**: 10 figures/tables
- **Review status**: Fully peer-reviewed, same standards as main figures
- **Referenced as**: "Extended Data Fig. 1", "Extended Data Table 1"
- **Use for**: Important results that support main findings but don't fit in the 6-item main limit
- **Examples**: Validation experiments, additional computational analyses, detailed breakdowns

### Supplementary Information (Not Fully Peer-Reviewed)
- **Maximum**: No limit
- **Review status**: Checked for completeness, NOT fully peer-reviewed
- **Referenced as**: "Supplementary Fig. 1", "Supplementary Table 1"
- **Use for**: Raw data, large tables, extended methods, additional analyses
- **Examples**: Full gene lists, complete parameter tables, raw quality metrics

### Decision Guide

| Content | Where It Goes |
|---|---|
| Key validation of a main result | Extended Data |
| Supporting analysis referenced in Discussion | Extended Data |
| Complete list of 127 genes (too large for main) | Supplementary Table |
| Full hyperparameter search results | Supplementary Table |
| Additional quality control metrics | Supplementary Figure |
| Raw benchmark comparisons | Supplementary Table |
| Detailed protocol steps | Supplementary Methods |

## Pre-Submission Figure Checklist

```
□ Total figures + tables ≤ 6
□ Extended Data items ≤ 10
□ All panels labeled (a, b, c) in lowercase bold
□ Font is Arial/Helvetica throughout
□ Font size ≥ 5 pt at final size
□ Resolution ≥ 300 DPI (≥ 600 for line art)
□ Colors are color-blind safe
□ All axes labeled with units
□ Scale bars included for microscopy
□ Legends ≤ 300 words each
□ Legends include n, test name, error bar definition
□ Every panel referenced in main text
□ Figures saved in TIFF/PDF/EPS format
□ No chart junk (3D, shadows, gradients)
□ Consistent style across all figures
□ Printed at column width — all labels readable
□ Extended Data figures meet same quality standards as main figures
```

