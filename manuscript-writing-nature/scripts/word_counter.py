"""
Word Counter for Scientific Manuscripts
=========================================
Counts words per section in LaTeX or DOCX manuscripts, comparing against
journal-specific limits.

Supports: LaTeX (.tex) and DOCX (.docx) files.

Usage:
    python word_counter.py manuscript/main.tex --profile nature
    python word_counter.py manuscript/paper.docx --profile nature
"""

import os
import re
import sys
import io
import argparse

# Fix Windows cp1252 encoding for Unicode output
if sys.stdout.encoding and sys.stdout.encoding.lower() != 'utf-8':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

os.chdir(os.path.dirname(os.path.abspath(__file__)))

# Import shared functions from quality_check
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from quality_check import (
    strip_latex_commands, extract_sections_latex,
    extract_sections_docx, count_words, PROFILES
)


def count_by_section(filepath, profile_name='nature'):
    """Count words in each section of a manuscript."""
    profile = PROFILES.get(profile_name)
    if not profile:
        print(f"Error: Unknown profile '{profile_name}'")
        sys.exit(1)

    is_docx = filepath.lower().endswith('.docx')
    if is_docx:
        sections = extract_sections_docx(filepath)
    else:
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
        sections = extract_sections_latex(content)

    print("=" * 60)
    print(f"WORD COUNT REPORT — {profile['name']} Profile")
    print(f"File: {filepath}")
    print("=" * 60)

    total_main = 0
    total_all = 0

    # Title (character count)
    title = sections.get('title', '')
    title_chars = len(title)
    title_limit = profile['title_max_chars']
    status = '✅' if title_chars <= title_limit else '❌'
    print(f"\n  {status} Title: {title_chars} chars (limit: {title_limit})")
    if title:
        print(f"       \"{title[:80]}{'...' if len(title) > 80 else ''}\"")

    # Summary/Abstract
    summary = sections.get('summary', sections.get('abstract', ''))
    summary_words = count_words(summary)
    summary_limit = profile['summary_max_words']
    status = '✅' if summary_words <= summary_limit else '❌'
    print(f"  {status} Summary: {summary_words} words (limit: {summary_limit})")
    total_all += summary_words

    # Main body sections
    print(f"\n  Main Body (limit: {profile['main_body_max_words']} words):")
    main_sections = ['introduction', 'results', 'discussion']
    for sec in main_sections:
        sec_text = sections.get(sec, '')
        sec_words = count_words(sec_text)
        total_main += sec_words
        total_all += sec_words
        print(f"    • {sec.title():15s}: {sec_words:5d} words")

    # If no individual sections found, use main_body
    if total_main == 0 and 'main_body' in sections:
        total_main = count_words(sections['main_body'])
        total_all += total_main
        print(f"    • Main body:       {total_main:5d} words (sections not detected)")

    status = '✅' if total_main <= profile['main_body_max_words'] else '❌'
    remaining = profile['main_body_max_words'] - total_main
    print(f"  {status} Main body total: {total_main} words", end='')
    if remaining >= 0:
        print(f" ({remaining} remaining)")
    else:
        print(f" ({-remaining} OVER LIMIT)")

    # Methods (no limit)
    methods = sections.get('methods', sections.get('online methods', ''))
    methods_words = count_words(methods)
    total_all += methods_words
    print(f"\n  ℹ️  Methods: {methods_words} words (no limit)")

    # Other sections
    print(f"\n  Other sections:")
    other_sections = ['acknowledgements', 'data_availability', 'code_availability',
                      'author_contributions', 'competing_interests']
    for sec in other_sections:
        for key, text in sections.items():
            if sec.replace('_', ' ') in key.replace('_', ' '):
                sec_words = count_words(text)
                total_all += sec_words
                print(f"    • {sec.replace('_', ' ').title():25s}: {sec_words:5d} words")
                break

    print(f"\n{'─' * 60}")
    print(f"  Total (all sections): {total_all} words")
    print(f"  Main body (counted): {total_main} / {profile['main_body_max_words']} words")
    print(f"{'─' * 60}")


def main():
    parser = argparse.ArgumentParser(
        description='Count words per section in scientific manuscripts'
    )
    parser.add_argument('manuscript', help='Path to manuscript file (.tex or .docx)')
    parser.add_argument('--profile', '-p', default='nature',
                        choices=list(PROFILES.keys()),
                        help='Journal profile (default: nature)')

    args = parser.parse_args()

    if not os.path.exists(args.manuscript):
        print(f"Error: File not found: {args.manuscript}")
        sys.exit(1)

    count_by_section(args.manuscript, args.profile)


if __name__ == '__main__':
    main()
