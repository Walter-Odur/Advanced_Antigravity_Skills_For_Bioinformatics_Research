"""
Citation Validator for Scientific Manuscripts
==============================================
Validates BibTeX references against the manuscript to ensure:
- All cited references exist in the .bib file
- All .bib entries are actually cited
- Reference count within journal limits
- Basic BibTeX entry completeness

Usage:
    python citation_validator.py references.bib main.tex --profile nature
    python citation_validator.py references.bib main.tex --output report.json
"""

import os
import re
import sys
import io
import json
import argparse

# Fix Windows cp1252 encoding for Unicode output
if sys.stdout.encoding and sys.stdout.encoding.lower() != 'utf-8':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

os.chdir(os.path.dirname(os.path.abspath(__file__)))


def parse_bib_file(bib_path):
    """Parse a BibTeX file and extract entries."""
    with open(bib_path, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()

    entries = {}
    # Match @type{key, ... }
    pattern = re.compile(
        r'@(\w+)\s*\{\s*([^,\s]+)\s*,\s*(.*?)\n\}',
        re.DOTALL
    )

    for match in pattern.finditer(content):
        entry_type = match.group(1).lower()
        key = match.group(2).strip()
        body = match.group(3)

        # Extract fields
        fields = {}
        field_pattern = re.compile(r'(\w+)\s*=\s*\{(.*?)\}', re.DOTALL)
        for field_match in field_pattern.finditer(body):
            field_name = field_match.group(1).lower()
            field_value = field_match.group(2).strip()
            fields[field_name] = field_value

        entries[key] = {
            'type': entry_type,
            'fields': fields,
        }

    return entries


def extract_citations_tex(tex_path):
    """Extract all citation keys used in a LaTeX file."""
    with open(tex_path, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()

    # Remove comments
    content = re.sub(r'%.*$', '', content, flags=re.MULTILINE)

    # Find all \cite variants
    cite_keys = set()
    for match in re.finditer(r'\\(?:cite|citep|citet|citealt|citealp|citeauthor|citeyear)\s*(?:\[[^\]]*\])?\s*\{([^}]+)\}', content):
        for key in match.group(1).split(','):
            cite_keys.add(key.strip())

    return cite_keys


def validate_entry_completeness(key, entry):
    """Check if a BibTeX entry has required fields."""
    issues = []
    entry_type = entry['type']
    fields = entry['fields']

    # Required fields by type
    required = {
        'article': ['author', 'title', 'journal', 'year'],
        'book': ['author', 'title', 'publisher', 'year'],
        'inproceedings': ['author', 'title', 'booktitle', 'year'],
        'incollection': ['author', 'title', 'booktitle', 'publisher', 'year'],
        'phdthesis': ['author', 'title', 'school', 'year'],
        'misc': ['author', 'title', 'year'],
        'techreport': ['author', 'title', 'institution', 'year'],
    }

    req_fields = required.get(entry_type, ['author', 'title', 'year'])
    for field in req_fields:
        if field not in fields or not fields[field].strip():
            issues.append(f"Missing required field: {field}")

    # Check for common issues
    if 'year' in fields:
        year = fields['year'].strip()
        if not re.match(r'^\d{4}$', year):
            issues.append(f"Invalid year format: '{year}'")

    if 'pages' in fields:
        pages = fields['pages']
        if '-' in pages and '--' not in pages:
            issues.append(f"Pages should use en-dash (--): '{pages}'")

    return issues


def run_validation(bib_path, tex_path, profile_name='nature'):
    """Run citation validation."""
    from quality_check import PROFILES
    profile = PROFILES.get(profile_name, PROFILES['nature'])

    print("=" * 60)
    print(f"CITATION VALIDATION — {profile['name']} Profile")
    print(f"BibTeX: {bib_path}")
    print(f"Manuscript: {tex_path}")
    print("=" * 60)

    # Parse files
    bib_entries = parse_bib_file(bib_path)
    cited_keys = extract_citations_tex(tex_path)

    print(f"\n  BibTeX entries: {len(bib_entries)}")
    print(f"  Citations in text: {len(cited_keys)}")
    print(f"  Reference limit: {profile['max_references']}")

    results = {
        'bib_entries': len(bib_entries),
        'cited_keys': len(cited_keys),
        'max_references': profile['max_references'],
        'issues': [],
    }

    # Check 1: Reference count
    print(f"\n{'─' * 40}")
    print("Check 1: Reference Count")
    if len(cited_keys) > profile['max_references']:
        print(f"  ❌ {len(cited_keys)} citations (max {profile['max_references']})")
        results['issues'].append({
            'type': 'over_limit',
            'severity': 'HIGH',
            'message': f"Too many references: {len(cited_keys)} (max {profile['max_references']})",
        })
    else:
        print(f"  ✅ {len(cited_keys)} citations (max {profile['max_references']})")

    # Check 2: Missing references (cited but not in .bib)
    print(f"\n{'─' * 40}")
    print("Check 2: Missing References")
    missing = cited_keys - set(bib_entries.keys())
    if missing:
        print(f"  ❌ {len(missing)} citations not found in .bib file:")
        for key in sorted(missing):
            print(f"     • {key}")
            results['issues'].append({
                'type': 'missing_bib',
                'severity': 'HIGH',
                'key': key,
                'message': f"Citation '{key}' not found in .bib file",
            })
    else:
        print(f"  ✅ All cited references found in .bib file")

    # Check 3: Uncited entries (in .bib but never cited)
    print(f"\n{'─' * 40}")
    print("Check 3: Uncited Entries")
    uncited = set(bib_entries.keys()) - cited_keys
    if uncited:
        print(f"  ⚠️  {len(uncited)} .bib entries never cited:")
        for key in sorted(uncited)[:10]:
            print(f"     • {key}")
        if len(uncited) > 10:
            print(f"     ... and {len(uncited) - 10} more")
        results['issues'].append({
            'type': 'uncited',
            'severity': 'LOW',
            'message': f"{len(uncited)} .bib entries are never cited",
            'keys': sorted(uncited),
        })
    else:
        print(f"  ✅ All .bib entries are cited")

    # Check 4: Entry completeness
    print(f"\n{'─' * 40}")
    print("Check 4: Entry Completeness")
    incomplete = 0
    for key, entry in bib_entries.items():
        entry_issues = validate_entry_completeness(key, entry)
        if entry_issues:
            incomplete += 1
            if incomplete <= 5:
                print(f"  ⚠️  {key}:")
                for issue in entry_issues:
                    print(f"     • {issue}")
            results['issues'].append({
                'type': 'incomplete',
                'severity': 'MEDIUM',
                'key': key,
                'details': entry_issues,
            })
    if incomplete > 5:
        print(f"  ... and {incomplete - 5} more entries with issues")
    if incomplete == 0:
        print(f"  ✅ All entries have required fields")

    # Summary
    high = sum(1 for i in results['issues'] if i.get('severity') == 'HIGH')
    med = sum(1 for i in results['issues'] if i.get('severity') == 'MEDIUM')
    low = sum(1 for i in results['issues'] if i.get('severity') == 'LOW')

    print(f"\n{'=' * 60}")
    print(f"SUMMARY: {len(results['issues'])} issues found")
    print(f"  High: {high}  Medium: {med}  Low: {low}")
    if high == 0:
        print(f"  ✅ No critical citation issues")
    else:
        print(f"  ❌ {high} critical issues must be fixed")
    print(f"{'=' * 60}")

    return results


def main():
    parser = argparse.ArgumentParser(
        description='Validate citations in scientific manuscripts'
    )
    parser.add_argument('bib', help='Path to BibTeX file (.bib)')
    parser.add_argument('tex', help='Path to manuscript file (.tex)')
    parser.add_argument('--profile', '-p', default='nature',
                        help='Journal profile (default: nature)')
    parser.add_argument('--output', '-o', default=None,
                        help='Save report as JSON file')

    args = parser.parse_args()

    for f in [args.bib, args.tex]:
        if not os.path.exists(f):
            print(f"Error: File not found: {f}")
            sys.exit(1)

    results = run_validation(args.bib, args.tex, args.profile)

    if args.output:
        with open(args.output, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"\nReport saved to: {args.output}")


if __name__ == '__main__':
    main()
