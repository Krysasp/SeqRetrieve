#!/usr/bin/env python3
"""
Generate an HTML report summarizing sequence retrieval, filtering, and classification.
"""

import sys
import os
import argparse
import csv
from typing import Dict, List
import re

def parse_year(date_str: str) -> str:
    if not date_str:
        return 'Unknown'
    match = re.search(r'(\d{4})', date_str)
    if match:
        return match.group(1)
    return 'Unknown'

def load_metadata(tsv_file: str) -> List[Dict]:
    """Load metadata TSV as list of dicts."""
    rows = []
    try:
        with open(tsv_file, 'r', newline='') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                rows.append(row)
    except Exception as e:
        print(f"Warning: Could not read metadata: {e}")
    return rows

def compute_classification(meta_rows: List[Dict]) -> Dict[str, Dict[str, int]]:
    """Compute counts for each classification field."""
    classification = {}
    fields = ['species', 'genus', 'organism', 'taxid', 'country', 'genotype', 'host']
    for field in fields:
        counts = {}
        for row in meta_rows:
            val = row.get(field, '').strip()
            if not val:
                val = 'Unknown'
            counts[val] = counts.get(val, 0) + 1
        # Sort by count descending
        sorted_counts = sorted(counts.items(), key=lambda x: x[1], reverse=True)
        classification[field] = dict(sorted_counts)
    # Year from collection_date
    year_counts = {}
    for row in meta_rows:
        date_str = row.get('collection_date', '')
        year = parse_year(date_str)
        year_counts[year] = year_counts.get(year, 0) + 1
    sorted_years = sorted(year_counts.items(), key=lambda x: x[0])
    classification['year'] = dict(sorted_years)
    return classification

def generate_html_report(
    output_html: str,
    fasta_file: str = None,
    metadata_file: str = None,
    filter_report: str = None,
    filter_stats: Dict = None,
    classification: Dict = None,
):
    """Generate HTML report."""
    # Read filter report if provided
    filtered_details = []
    if filter_report and os.path.exists(filter_report):
        try:
            with open(filter_report, 'r', newline='') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    filtered_details.append(row)
        except Exception:
            pass

    # Generate HTML
    from datetime import datetime
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    html = f"""<!DOCTYPE html>
<html>
<head>
    <title>RetrieveSeq Report</title>
    <meta charset="utf-8">
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; color: #333; }}
        h1 {{ color: #2c3e50; }}
        h2 {{ color: #34495e; border-bottom: 1px solid #ccc; padding-bottom: 5px; }}
        table {{ border-collapse: collapse; width: 100%; margin-bottom: 20px; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #f2f2f2; }}
        tr:nth-child(even) {{ background-color: #f9f9f9; }}
        .summary {{ background-color: #e8f4f8; padding: 15px; border-radius: 5px; margin-bottom: 20px; }}
        .filter-section {{ background-color: #fff3cd; padding: 15px; border-radius: 5px; margin-bottom: 20px; }}
        .class-section {{ background-color: #d4edda; padding: 15px; border-radius: 5px; margin-bottom: 20px; }}
        .count {{ font-weight: bold; color: #2980b9; }}
        .bar {{ background-color: #3498db; height: 20px; margin: 2px 0; }}
    </style>
</head>
<body>
    <h1>RetrieveSeq Analysis Report</h1>
    <p>Generated on {timestamp}</p>

    <div class="summary">
        <h2>Summary</h2>
        <p>FASTA file: {fasta_file or 'N/A'}</p>
        <p>Metadata file: {metadata_file or 'N/A'}</p>
"""
    # Add filter stats
    if filter_stats:
        html += f"""
        <p>Total sequences: <span class="count">{filter_stats.get('total', 'N/A')}</span></p>
        <p>Passed filters: <span class="count">{filter_stats.get('passed', 'N/A')}</span></p>
        <p>Filtered out: <span class="count">{filter_stats.get('filtered', 'N/A')}</span></p>
        <p>Longest sequence: <span class="count">{filter_stats.get('longest', 'N/A'):,}</span> bp</p>
"""
    html += """
    </div>
"""

    # Filter details
    if filtered_details:
        html += """
    <div class="filter-section">
        <h2>Filtering Details</h2>
        <p>Sequences filtered out:</p>
        <table>
            <tr><th>Accession</th><th>Length</th><th>N count</th><th>Total ambig</th><th>Ambig %</th><th>Reasons</th></tr>
"""
        for row in filtered_details[:100]:  # limit to first 100
            html += f"            <tr><td>{row.get('accession', '')}</td><td>{row.get('length', '')}</td><td>{row.get('n_count', '')}</td><td>{row.get('total_ambig', '')}</td><td>{row.get('ambig_pct', '')}</td><td>{row.get('filter_reasons', '')}</td></tr>\n"
        if len(filtered_details) > 100:
            html += f"            <tr><td colspan=6>... and {len(filtered_details)-100} more</td></tr>\n"
        html += """
        </table>
    </div>
"""

    # Classification
    if classification:
        html += """
    <div class="class-section">
        <h2>Classification by Field</h2>
"""
        for field, counts in classification.items():
            html += f"""
        <h3>{field.title()} (top 20)</h3>
        <table>
            <tr><th>{field.title()}</th><th>Count</th><th>Percentage</th></tr>
"""
            total = sum(counts.values())
            for val, cnt in list(counts.items())[:20]:
                pct = cnt / total * 100 if total > 0 else 0
                html += f"            <tr><td>{val}</td><td>{cnt:,}</td><td>{pct:.1f}%</td></tr>\n"
            if len(counts) > 20:
                html += f"            <tr><td colspan=3>... and {len(counts)-20} more categories</td></tr>\n"
            html += """
        </table>
"""
        html += """
    </div>
"""

    html += """
</body>
</html>
"""
    # Write file
    with open(output_html, 'w') as f:
        f.write(html)
    print(f"HTML report written to {output_html}")

def main():
    parser = argparse.ArgumentParser(
        description='Generate HTML report for RetrieveSeq pipeline.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python generate_report.py -o report.html -f seqs.fasta -m meta.tsv --filter-report filter_report.tsv
        """
    )
    parser.add_argument('-o', '--output', required=True, help='Output HTML file')
    parser.add_argument('-f', '--fasta', default=None, help='Input FASTA file (for stats)')
    parser.add_argument('-m', '--metadata', default=None, help='Metadata TSV file')
    parser.add_argument('--filter-report', default=None, help='Filter report TSV file')
    parser.add_argument('--filter-stats', default=None, help='JSON file with filter stats (optional)')
    args = parser.parse_args()

    # Load metadata
    meta_rows = []
    if args.metadata and os.path.exists(args.metadata):
        meta_rows = load_metadata(args.metadata)
        print(f"Loaded {len(meta_rows):,} metadata records")

    # Compute classification
    classification = None
    if meta_rows:
        classification = compute_classification(meta_rows)

    # Load filter stats if JSON provided
    filter_stats = None
    if args.filter_stats and os.path.exists(args.filter_stats):
        import json
        with open(args.filter_stats, 'r') as f:
            filter_stats = json.load(f)

    # Generate report
    generate_html_report(
        output_html=args.output,
        fasta_file=args.fasta,
        metadata_file=args.metadata,
        filter_report=args.filter_report,
        filter_stats=filter_stats,
        classification=classification,
    )

if __name__ == '__main__':
    main()
