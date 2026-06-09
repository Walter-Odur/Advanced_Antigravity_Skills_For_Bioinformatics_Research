"""
Automated 10-Pass Quality Audit Engine for Scientific Manuscripts
==================================================================
Performs automated checks on LaTeX manuscripts against journal-specific
profiles. Generates a structured report with pass/fail status, severity
ratings, and specific fix suggestions.

Supports: LaTeX (.tex) and DOCX (.docx) manuscripts.

Usage:
    python quality_check.py manuscript/main.tex --profile nature
    python quality_check.py manuscript/paper.docx --profile nature
    python quality_check.py manuscript/main.tex --profile nature --output report.json
"""

import os
import re
import sys
import io
import json
import argparse
from pathlib import Path
from collections import Counter

# Fix Windows cp1252 encoding for Unicode output
if sys.stdout.encoding and sys.stdout.encoding.lower() != 'utf-8':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

os.chdir(os.path.dirname(os.path.abspath(__file__)))


# ── Journal Profiles ──────────────────────────────────────────────
PROFILES = {
    'nature': {
        'name': 'Nature',
        'summary_max_words': 200,
        'main_body_max_words': 3000,
        'legend_max_words': 300,
        'title_max_chars': 90,
        'max_references': 50,
        'max_display_items': 6,
        'max_extended_data': 10,
        'required_sections': [
            'data_availability', 'code_availability',
            'author_contributions', 'competing_interests',
        ],
        'summary_label': 'Summary',
        'uses_section_headers': False,
    },
}


# ── LaTeX Text Extraction ────────────────────────────────────────
def strip_latex_commands(text):
    """Remove LaTeX commands, leaving readable text for word counting."""
    # Remove comments
    text = re.sub(r'%.*$', '', text, flags=re.MULTILINE)
    # Remove common environments we don't count
    for env in ['equation', 'equation*', 'align', 'align*', 'figure',
                'table', 'tikzpicture', 'lstlisting']:
        text = re.sub(
            rf'\\begin\{{{env}\}}.*?\\end\{{{env}\}}',
            '', text, flags=re.DOTALL
        )
    # Remove \cite, \ref, \label commands
    text = re.sub(r'\\(?:cite|ref|label|eqref|autoref|nameref)\{[^}]*\}', '', text)
    # Remove \textbf, \textit, \emph — keep content
    text = re.sub(r'\\(?:textbf|textit|emph|underline)\{([^}]*)\}', r'\1', text)
    # Remove remaining commands
    text = re.sub(r'\\[a-zA-Z]+\*?(?:\[[^\]]*\])?(?:\{[^}]*\})?', '', text)
    # Remove braces
    text = re.sub(r'[{}]', '', text)
    # Clean whitespace
    text = re.sub(r'\s+', ' ', text).strip()
    return text


def extract_sections_latex(content):
    """Extract logical sections from LaTeX content."""
    sections = {}

    # Extract title
    title_match = re.search(r'\\title\{([^}]+)\}', content)
    if title_match:
        sections['title'] = strip_latex_commands(title_match.group(1))

    # Extract abstract/summary
    abstract_match = re.search(
        r'\\begin\{abstract\}(.*?)\\end\{abstract\}',
        content, re.DOTALL
    )
    if not abstract_match:
        abstract_match = re.search(
            r'\\abstract\{(.*?)\}',
            content, re.DOTALL
        )
    if abstract_match:
        sections['summary'] = strip_latex_commands(abstract_match.group(1))

    # Extract sections by \section commands
    section_pattern = re.compile(
        r'\\section\*?\{([^}]+)\}(.*?)(?=\\section\*?\{|\\begin\{thebibliography\}|\\bibliography\{|$)',
        re.DOTALL
    )
    for match in section_pattern.finditer(content):
        name = match.group(1).strip().lower()
        text = strip_latex_commands(match.group(2))
        sections[name] = text

    # Full body text (everything between abstract and methods/bibliography)
    body_match = re.search(
        r'\\end\{abstract\}(.*?)(?=\\section\*?\{[Mm]ethod|\\begin\{thebibliography\}|\\bibliography\{)',
        content, re.DOTALL
    )
    if body_match:
        sections['main_body'] = strip_latex_commands(body_match.group(1))

    return sections


def extract_sections_docx(filepath):
    """Extract sections from a DOCX file."""
    sections = {}
    try:
        import zipfile
        from xml.etree import ElementTree as ET

        ns = {'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}

        with zipfile.ZipFile(filepath) as z:
            with z.open('word/document.xml') as f:
                tree = ET.parse(f)
                root = tree.getroot()

        paragraphs = root.findall('.//w:p', ns)
        full_text = []
        current_section = 'unknown'
        section_texts = {}

        for para in paragraphs:
            # Get paragraph style
            style_el = para.find('.//w:pStyle', ns)
            style = style_el.get(f'{{{ns["w"]}}}val', '') if style_el is not None else ''

            # Get text content
            texts = para.findall('.//w:t', ns)
            para_text = ''.join(t.text or '' for t in texts)

            if not para_text.strip():
                continue

            # Detect section headers by style
            if 'Heading' in style or 'Title' in style:
                heading_lower = para_text.strip().lower()
                if 'abstract' in heading_lower or 'summary' in heading_lower:
                    current_section = 'summary'
                elif 'introduction' in heading_lower:
                    current_section = 'introduction'
                elif 'result' in heading_lower:
                    current_section = 'results'
                elif 'discussion' in heading_lower:
                    current_section = 'discussion'
                elif 'method' in heading_lower:
                    current_section = 'methods'
                elif 'data availability' in heading_lower:
                    current_section = 'data_availability'
                elif 'code availability' in heading_lower:
                    current_section = 'code_availability'
                elif 'author contribution' in heading_lower:
                    current_section = 'author_contributions'
                elif 'competing' in heading_lower or 'conflict' in heading_lower:
                    current_section = 'competing_interests'
                elif 'acknowledge' in heading_lower:
                    current_section = 'acknowledgements'
                elif 'reference' in heading_lower:
                    current_section = 'references'
                else:
                    current_section = heading_lower
                if 'Title' in style:
                    sections['title'] = para_text.strip()
                continue

            if current_section not in section_texts:
                section_texts[current_section] = []
            section_texts[current_section].append(para_text.strip())
            full_text.append(para_text.strip())

        for sec, texts in section_texts.items():
            sections[sec] = ' '.join(texts)

        # Build main_body from intro + results + discussion
        main_parts = []
        for key in ['introduction', 'results', 'discussion']:
            if key in sections:
                main_parts.append(sections[key])
        if main_parts:
            sections['main_body'] = ' '.join(main_parts)

    except Exception as e:
        print(f"  Warning: DOCX parsing error: {e}")

    return sections


def count_words(text):
    """Count words in text."""
    if not text:
        return 0
    words = text.split()
    return len(words)


def count_references_latex(content):
    """Count references in LaTeX."""
    # Count \bibitem entries
    bibitems = re.findall(r'\\bibitem', content)
    if bibitems:
        return len(bibitems)

    # Count entries in .bib file reference
    bib_match = re.search(r'\\bibliography\{([^}]+)\}', content)
    if bib_match:
        bib_file = bib_match.group(1)
        if not bib_file.endswith('.bib'):
            bib_file += '.bib'
        # Try to find and count entries in .bib file
        for search_dir in ['.', '..', '../manuscript']:
            bib_path = os.path.join(search_dir, bib_file)
            if os.path.exists(bib_path):
                with open(bib_path, 'r', encoding='utf-8', errors='ignore') as f:
                    bib_content = f.read()
                return len(re.findall(r'@\w+\{', bib_content))

    # Count unique citation keys used
    cite_keys = set()
    for match in re.finditer(r'\\cite[tp]?\{([^}]+)\}', content):
        for key in match.group(1).split(','):
            cite_keys.add(key.strip())
    return len(cite_keys)


def count_display_items_latex(content):
    """Count figures and tables in LaTeX."""
    figures = len(re.findall(r'\\begin\{figure\*?\}', content))
    tables = len(re.findall(r'\\begin\{table\*?\}', content))
    return figures, tables


def find_figure_legends_latex(content):
    """Extract figure legends and count their words."""
    legends = []
    for match in re.finditer(
        r'\\begin\{figure\*?\}.*?\\caption\{(.*?)\}.*?\\end\{figure\*?\}',
        content, re.DOTALL
    ):
        legend_text = strip_latex_commands(match.group(1))
        legends.append(legend_text)
    return legends


def check_acronyms(text):
    """Find undefined acronyms (uppercase sequences ≥2 chars)."""
    acronyms = re.findall(r'\b[A-Z]{2,}\b', text)
    acronym_counts = Counter(acronyms)
    # Check if each acronym is defined (word followed by parenthetical acronym)
    undefined = []
    for acronym, count in acronym_counts.items():
        if acronym in ('AND', 'THE', 'FOR', 'NOT', 'BUT', 'OR', 'NOR', 'YET',
                       'DNA', 'RNA', 'PCR', 'USA', 'UK', 'WHO', 'HIV', 'AIDS',
                       'TB', 'COVID', 'SARS', 'UMAP', 'PCA', 'ROC', 'AUC',
                       'CI', 'SD', 'IQR', 'FDR', 'CPU', 'GPU', 'RAM', 'DOI',
                       'FAIR', 'ORCID', 'PDF', 'TIFF', 'EPS', 'CSV', 'JSON',
                       'HTML', 'HTTP', 'HTTPS', 'API', 'URL', 'GEO', 'SRA'):
            continue
        # Check for definition pattern: "full name (ACRONYM)"
        pattern = rf'\([^)]*{re.escape(acronym)}[^)]*\)'
        if not re.search(pattern, text):
            undefined.append(acronym)
    return undefined


def check_sentence_length(text):
    """Find sentences exceeding 40 words."""
    sentences = re.split(r'[.!?]+', text)
    long_sentences = []
    for sent in sentences:
        words = sent.strip().split()
        if len(words) > 40:
            preview = ' '.join(words[:15]) + '...'
            long_sentences.append({
                'words': len(words),
                'preview': preview,
            })
    return long_sentences


def check_passive_voice(text):
    """Detect passive voice constructions."""
    passive_patterns = [
        r'\bwas\s+\w+ed\b', r'\bwere\s+\w+ed\b',
        r'\bis\s+\w+ed\b', r'\bare\s+\w+ed\b',
        r'\bbeen\s+\w+ed\b', r'\bwas\s+\w+en\b',
        r'\bwere\s+\w+en\b',
    ]
    count = 0
    for pattern in passive_patterns:
        count += len(re.findall(pattern, text, re.IGNORECASE))
    return count


def check_todo_comments(content):
    """Find TODO/FIXME/XXX comments."""
    todos = re.findall(r'(?:TODO|FIXME|XXX|HACK|PLACEHOLDER)[:\s].*$',
                       content, re.MULTILINE | re.IGNORECASE)
    return todos


# ── Main Quality Check Engine ────────────────────────────────────
def run_quality_check(filepath, profile_name='nature'):
    """Run the full 10-pass quality check."""
    profile = PROFILES.get(profile_name)
    if not profile:
        print(f"Error: Unknown profile '{profile_name}'")
        print(f"Available profiles: {', '.join(PROFILES.keys())}")
        sys.exit(1)

    print("=" * 70)
    print(f"MANUSCRIPT QUALITY AUDIT — {profile['name']} Profile")
    print(f"File: {filepath}")
    print("=" * 70)

    # Read file
    is_docx = filepath.lower().endswith('.docx')
    if is_docx:
        content = ''  # DOCX doesn't have raw text content in the same way
        sections = extract_sections_docx(filepath)
    else:
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
        sections = extract_sections_latex(content)

    results = {
        'profile': profile_name,
        'file': filepath,
        'passes': {},
        'summary': {'passed': 0, 'failed': 0, 'warnings': 0},
    }

    # ── Pass 1: Structure & Completeness ──
    print(f"\n{'─' * 50}")
    print("PASS 1: Structure & Completeness")
    print(f"{'─' * 50}")
    pass1 = {'status': 'PASS', 'issues': []}

    # Check required sections
    for req in profile['required_sections']:
        req_found = False
        for key in sections:
            if req.replace('_', ' ') in key.replace('_', ' '):
                req_found = True
                break
        # Also check in raw content for LaTeX
        if not req_found and not is_docx:
            patterns = {
                'data_availability': r'[Dd]ata\s+[Aa]vailability',
                'code_availability': r'[Cc]ode\s+[Aa]vailability',
                'author_contributions': r'[Aa]uthor\s+[Cc]ontribution',
                'competing_interests': r'[Cc]ompeting\s+[Ii]nterest|[Cc]onflict\s+of\s+[Ii]nterest',
            }
            if req in patterns:
                if re.search(patterns[req], content):
                    req_found = True
        status = '✅' if req_found else '❌'
        if not req_found:
            pass1['issues'].append({
                'severity': 'HIGH',
                'message': f"Missing required section: {req.replace('_', ' ').title()}",
            })
            pass1['status'] = 'FAIL'
        print(f"  {status} {req.replace('_', ' ').title()}")

    # Check for TODO comments
    if not is_docx:
        todos = check_todo_comments(content)
        if todos:
            pass1['issues'].append({
                'severity': 'MEDIUM',
                'message': f"Found {len(todos)} TODO/FIXME comments",
                'details': todos[:5],
            })
            print(f"  ⚠️  {len(todos)} TODO/FIXME comments found")
        else:
            print(f"  ✅ No TODO/FIXME comments")

    results['passes']['1_structure'] = pass1
    print(f"  → Pass 1: {pass1['status']}")

    # ── Pass 2: Word Count Compliance ──
    print(f"\n{'─' * 50}")
    print("PASS 2: Word Count Compliance")
    print(f"{'─' * 50}")
    pass2 = {'status': 'PASS', 'issues': [], 'counts': {}}

    # Title
    title = sections.get('title', '')
    title_chars = len(title)
    pass2['counts']['title_chars'] = title_chars
    if title_chars > profile['title_max_chars']:
        pass2['issues'].append({
            'severity': 'HIGH',
            'message': f"Title too long: {title_chars} chars (max {profile['title_max_chars']})",
        })
        pass2['status'] = 'FAIL'
        print(f"  ❌ Title: {title_chars} chars (max {profile['title_max_chars']})")
    elif title:
        print(f"  ✅ Title: {title_chars} chars (max {profile['title_max_chars']})")

    # Summary/Abstract
    summary = sections.get('summary', sections.get('abstract', ''))
    summary_words = count_words(summary)
    pass2['counts']['summary_words'] = summary_words
    if summary_words > profile['summary_max_words']:
        pass2['issues'].append({
            'severity': 'HIGH',
            'message': f"Summary too long: {summary_words} words (max {profile['summary_max_words']})",
        })
        pass2['status'] = 'FAIL'
        print(f"  ❌ Summary: {summary_words} words (max {profile['summary_max_words']})")
    elif summary:
        print(f"  ✅ Summary: {summary_words} words (max {profile['summary_max_words']})")

    # Main body
    main_body = sections.get('main_body', '')
    body_words = count_words(main_body)
    pass2['counts']['main_body_words'] = body_words
    if body_words > profile['main_body_max_words']:
        over = body_words - profile['main_body_max_words']
        pass2['issues'].append({
            'severity': 'HIGH',
            'message': f"Main body too long: {body_words} words (max {profile['main_body_max_words']}, over by {over})",
        })
        pass2['status'] = 'FAIL'
        print(f"  ❌ Main body: {body_words} words (max {profile['main_body_max_words']}, +{over})")
    elif main_body:
        print(f"  ✅ Main body: {body_words} words (max {profile['main_body_max_words']})")

    # Figure legends
    if not is_docx:
        legends = find_figure_legends_latex(content)
        for i, legend in enumerate(legends, 1):
            legend_words = count_words(legend)
            if legend_words > profile['legend_max_words']:
                pass2['issues'].append({
                    'severity': 'MEDIUM',
                    'message': f"Fig. {i} legend too long: {legend_words} words (max {profile['legend_max_words']})",
                })
                pass2['status'] = 'FAIL'
                print(f"  ❌ Fig. {i} legend: {legend_words} words (max {profile['legend_max_words']})")
            else:
                print(f"  ✅ Fig. {i} legend: {legend_words} words")

    results['passes']['2_word_count'] = pass2
    print(f"  → Pass 2: {pass2['status']}")

    # ── Pass 5: Figures & Tables ──
    print(f"\n{'─' * 50}")
    print("PASS 5: Figures & Tables")
    print(f"{'─' * 50}")
    pass5 = {'status': 'PASS', 'issues': []}

    if not is_docx:
        n_figs, n_tabs = count_display_items_latex(content)
        total_display = n_figs + n_tabs
        pass5['counts'] = {'figures': n_figs, 'tables': n_tabs, 'total': total_display}

        if total_display > profile['max_display_items']:
            pass5['issues'].append({
                'severity': 'HIGH',
                'message': f"Too many display items: {total_display} (max {profile['max_display_items']})",
            })
            pass5['status'] = 'FAIL'
            print(f"  ❌ Display items: {total_display} (max {profile['max_display_items']})")
        else:
            print(f"  ✅ Display items: {total_display} (max {profile['max_display_items']})")

        # Check figure references in text
        fig_refs = set(re.findall(r'Fig(?:ure)?\.?\s*(\d+)', content))
        print(f"  ℹ️  Figures referenced in text: {sorted(fig_refs)}")

    results['passes']['5_figures'] = pass5
    print(f"  → Pass 5: {pass5['status']}")

    # ── Pass 7: Citation Integrity ──
    print(f"\n{'─' * 50}")
    print("PASS 7: Citation Integrity")
    print(f"{'─' * 50}")
    pass7 = {'status': 'PASS', 'issues': []}

    if not is_docx:
        n_refs = count_references_latex(content)
        pass7['counts'] = {'references': n_refs}

        if n_refs > profile['max_references']:
            pass7['issues'].append({
                'severity': 'HIGH',
                'message': f"Too many references: {n_refs} (max {profile['max_references']})",
            })
            pass7['status'] = 'FAIL'
            print(f"  ❌ References: {n_refs} (max {profile['max_references']})")
        else:
            print(f"  ✅ References: {n_refs} (max {profile['max_references']})")

    results['passes']['7_citations'] = pass7
    print(f"  → Pass 7: {pass7['status']}")

    # ── Pass 4: Accessibility & Jargon ──
    print(f"\n{'─' * 50}")
    print("PASS 4: Accessibility & Jargon")
    print(f"{'─' * 50}")
    pass4 = {'status': 'PASS', 'issues': []}

    # Check undefined acronyms
    full_text = ' '.join(sections.values())
    undefined = check_acronyms(full_text)
    if undefined:
        pass4['issues'].append({
            'severity': 'LOW',
            'message': f"Potentially undefined acronyms: {', '.join(undefined[:10])}",
        })
        print(f"  ⚠️  Potentially undefined: {', '.join(undefined[:10])}")
    else:
        print(f"  ✅ No undefined acronyms detected")

    # Check sentence length
    long_sents = check_sentence_length(full_text)
    if long_sents:
        pass4['issues'].append({
            'severity': 'LOW',
            'message': f"{len(long_sents)} sentences exceed 40 words",
            'details': long_sents[:3],
        })
        print(f"  ⚠️  {len(long_sents)} sentences exceed 40 words")
    else:
        print(f"  ✅ All sentences within length guidelines")

    # Check passive voice density
    passive_count = check_passive_voice(full_text)
    total_sentences = len(re.split(r'[.!?]+', full_text))
    if total_sentences > 0:
        passive_pct = 100 * passive_count / max(1, total_sentences)
        if passive_pct > 40:
            pass4['issues'].append({
                'severity': 'MEDIUM',
                'message': f"High passive voice usage: {passive_pct:.0f}% of sentences",
            })
            print(f"  ⚠️  Passive voice: {passive_pct:.0f}% (aim for <30%)")
        else:
            print(f"  ✅ Passive voice: {passive_pct:.0f}%")

    results['passes']['4_accessibility'] = pass4
    print(f"  → Pass 4: {pass4['status']}")

    # ── Summary ──
    total_issues = sum(len(p.get('issues', [])) for p in results['passes'].values())
    high_issues = sum(
        1 for p in results['passes'].values()
        for issue in p.get('issues', [])
        if issue.get('severity') == 'HIGH'
    )
    failed = sum(1 for p in results['passes'].values() if p['status'] == 'FAIL')
    passed = len(results['passes']) - failed

    results['summary'] = {
        'passed': passed,
        'failed': failed,
        'total_issues': total_issues,
        'high_severity': high_issues,
    }

    print(f"\n{'=' * 70}")
    print(f"AUDIT SUMMARY")
    print(f"{'=' * 70}")
    print(f"  Passes completed: {len(results['passes'])}")
    print(f"  Passed: {passed}")
    print(f"  Failed: {failed}")
    print(f"  Total issues: {total_issues}")
    print(f"  High severity: {high_issues}")

    if failed == 0:
        print(f"\n  ✅ ALL AUTOMATED CHECKS PASSED")
    else:
        print(f"\n  ❌ {failed} PASS(ES) FAILED — see details above")

    print(f"\n  Note: Passes 3, 6, 8, 9, 10 require manual review.")
    print(f"  See references/quality_checklist.md for the full protocol.")
    print(f"{'=' * 70}")

    return results


def main():
    parser = argparse.ArgumentParser(
        description='Automated quality audit for scientific manuscripts'
    )
    parser.add_argument('manuscript', help='Path to manuscript file (.tex or .docx)')
    parser.add_argument('--profile', '-p', default='nature',
                        choices=list(PROFILES.keys()),
                        help='Journal profile (default: nature)')
    parser.add_argument('--output', '-o', default=None,
                        help='Save report as JSON file')

    args = parser.parse_args()

    if not os.path.exists(args.manuscript):
        print(f"Error: File not found: {args.manuscript}")
        sys.exit(1)

    results = run_quality_check(args.manuscript, args.profile)

    if args.output:
        with open(args.output, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"\nReport saved to: {args.output}")


if __name__ == '__main__':
    main()
