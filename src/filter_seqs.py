#!/usr/bin/env python3
"""
Filter sequences based on length (percentage of longest), ambiguous base count/percentage.
Outputs filtered sequences to a subfolder, and passes to main output.
Generates statistics for reporting.
"""

import sys
import os
import argparse
import csv
from Bio import SeqIO
from typing import List, Dict, Tuple
import re

IUPAC_AMBIG = set('RYSWKMBDHV')

def compute_seq_stats(seq_record):
    """Compute statistics for a sequence."""
    seq_str = str(seq_record.seq).upper()
    length = len(seq_str)
    n_count = seq_str.count('N')
    # Count other IUPAC ambiguous codes (excluding ACGTN)
    other_ambig = sum(1 for c in seq_str if c in IUPAC_AMBIG)
    total_ambig = n_count + other_ambig
    ambig_pct = (total_ambig / length * 100) if length > 0 else 0.0
    n_pct = (n_count / length * 100) if length > 0 else 0.0
    return {
        'length': length,
        'n_count': n_count,
        'other_ambig': other_ambig,
        'total_ambig': total_ambig,
        'ambig_pct': ambig_pct,
        'n_pct': n_pct,
    }

def filter_sequences(
    fasta_file: str,
    output_dir: str,
    filtered_subdir: str,
    min_length_pct: float = None,
    max_ambig_count: int = None,
    max_ambig_pct: float = None,
) -> Tuple[List, List, Dict]:
    """
    Filter sequences.
    Returns: (passed_records, filtered_records_with_reasons, stats_dict)
    """
    print(f"Reading sequences from {fasta_file}...")
    records = list(SeqIO.parse(fasta_file, "fasta"))
    print(f"Total sequences: {len(records):,}")

    # Compute stats for all sequences
    seq_stats = []
    for rec in records:
        stats = compute_seq_stats(rec)
        seq_stats.append((rec, stats))

    # Determine longest length
    longest = max(stats['length'] for _, stats in seq_stats)
    print(f"Longest sequence length: {longest:,} bp")

    # Apply filters
    passed = []
    filtered = []
    filter_reasons = {'length': 0, 'ambig_count': 0, 'ambig_pct': 0}

    for rec, stats in seq_stats:
        fail_reasons = []
        # Length filter
        if min_length_pct is not None:
            min_len = longest * (min_length_pct / 100.0)
            if stats['length'] < min_len:
                fail_reasons.append(f"length < {min_len:.0f}")
                filter_reasons['length'] += 1
        # Ambig count filter
        if max_ambig_count is not None:
            if stats['total_ambig'] > max_ambig_count:
                fail_reasons.append(f"ambig_count > {max_ambig_count}")
                filter_reasons['ambig_count'] += 1
        # Ambig pct filter
        if max_ambig_pct is not None:
            if stats['ambig_pct'] > max_ambig_pct:
                fail_reasons.append(f"ambig_pct > {max_ambig_pct:.1f}%")
                filter_reasons['ambig_pct'] += 1

        if fail_reasons:
            filtered.append((rec, stats, fail_reasons))
        else:
            passed.append((rec, stats))

    print(f"\nFiltering results:")
    print(f"  Passed: {len(passed):,}")
    print(f"  Filtered out: {len(filtered):,}")
    if filter_reasons['length']:
        print(f"    - Failed length filter: {filter_reasons['length']:,}")
    if filter_reasons['ambig_count']:
        print(f"    - Failed ambig count filter: {filter_reasons['ambig_count']:,}")
    if filter_reasons['ambig_pct']:
        print(f"    - Failed ambig pct filter: {filter_reasons['ambig_pct']:,}")

    # Write passed sequences
    os.makedirs(output_dir, exist_ok=True)
    passed_file = os.path.join(output_dir, os.path.basename(fasta_file))
    with open(passed_file, 'w') as f:
        for rec, _ in passed:
            SeqIO.write(rec, f, "fasta")
    print(f"\nPassed sequences written to {passed_file}")

    # Write filtered sequences to subdir
    filtered_dir = os.path.join(output_dir, filtered_subdir)
    os.makedirs(filtered_dir, exist_ok=True)
    filtered_file = os.path.join(filtered_dir, os.path.basename(fasta_file))
    with open(filtered_file, 'w') as f:
        for rec, _, _ in filtered:
            SeqIO.write(rec, f, "fasta")
    print(f"Filtered sequences written to {filtered_file}")

    # Write filter report TSV
    report_file = os.path.join(output_dir, "filter_report.tsv")
    with open(report_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['accession', 'length', 'n_count', 'total_ambig', 'ambig_pct', 'filter_reasons'])
        for rec, stats, reasons in filtered:
            acc = rec.id.split()[0]
            writer.writerow([acc, stats['length'], stats['n_count'], stats['total_ambig'],
                           f"{stats['ambig_pct']:.2f}", "; ".join(reasons)])
    print(f"Filter report written to {report_file}")

    # Compute summary stats for report generation
    stats_summary = {
        'total': len(records),
        'passed': len(passed),
        'filtered': len(filtered),
        'longest': longest,
        'filter_reasons': filter_reasons,
        'passed_file': passed_file,
        'filtered_file': filtered_file,
        'report_file': report_file,
    }
    return passed, filtered, stats_summary

def main():
    parser = argparse.ArgumentParser(
        description='Filter sequences based on length and ambiguous base content.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Filter by length >= 80% of longest, max 5% ambiguous bases
  python filter_seqs.py -f seqs.fasta -o output/ --min-length-pct 80 --max-ambig-pct 5

  # Filter by max 100 ambiguous bases (N or other IUPAC)
  python filter_seqs.py -f seqs.fasta -o output/ --max-ambig-count 100

  # Combined filters
  python filter_seqs.py -f seqs.fasta -o output/ --min-length-pct 70 --max-ambig-count 50 --max-ambig-pct 2
        """
    )
    parser.add_argument('-f', '--fasta', required=True, help='Input FASTA file')
    parser.add_argument('-o', '--output-dir', required=True, help='Output directory for passed sequences')
    parser.add_argument('--filtered-subdir', default='filtered',
                        help='Subdirectory name for filtered sequences (default: filtered)')
    parser.add_argument('--min-length-pct', type=float, default=None,
                        help='Minimum length as percentage of longest sequence (e.g., 80.0)')
    parser.add_argument('--max-ambig-count', type=int, default=None,
                        help='Maximum count of ambiguous bases (N + other IUPAC codes)')
    parser.add_argument('--max-ambig-pct', type=float, default=None,
                        help='Maximum percentage of ambiguous bases')
    args = parser.parse_args()

    if not os.path.exists(args.fasta):
        print(f"Error: FASTA file {args.fasta} not found")
        sys.exit(1)

    passed, filtered, stats = filter_sequences(
        args.fasta,
        args.output_dir,
        args.filtered_subdir,
        min_length_pct=args.min_length_pct,
        max_ambig_count=args.max_ambig_count,
        max_ambig_pct=args.max_ambig_pct,
    )

    print(f"\nDone. Passed: {stats['passed']:,}, Filtered: {stats['filtered']:,}")

if __name__ == '__main__':
    main()
