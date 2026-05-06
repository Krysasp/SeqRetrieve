#!/usr/bin/env python3
"""
Sort sequences (FASTA) and metadata (TSV) by a specified field.
Supports: taxid, organism, species, genus, collection_date (year),
or any column present in the metadata TSV.
"""

import sys
import os
import argparse
import csv
from Bio import SeqIO
from typing import List, Dict, Optional

def parse_year(date_str: str) -> str:
    """Extract year from collection date string."""
    if not date_str:
        return '0000'
    # Handle formats like '2023', '2023-01', '2023-01-15', 'Jan-2023'
    import re
    match = re.search(r'(\d{4})', date_str)
    if match:
        return match.group(1)
    return '0000'

def load_metadata(tsv_file: str) -> Dict[str, dict]:
    """Load metadata TSV, keyed by accession (with and without version)."""
    meta = {}
    try:
        with open(tsv_file, 'r', newline='') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                acc = row.get('accession', '').strip()
                ver = row.get('version', '').strip()
                if acc:
                    meta[acc] = row
                if ver:
                    meta[ver] = row
                    # Also add without version
                    base = ver.split('.')[0]
                    if base not in meta:
                        meta[base] = row
    except Exception as e:
        print(f"Error reading metadata TSV: {e}")
        sys.exit(1)
    return meta

def load_sequences(fasta_file: str):
    """Load FASTA sequences, return list of (record, accession_base)."""
    seqs = []
    try:
        for record in SeqIO.parse(fasta_file, "fasta"):
            # Extract accession from record.id (e.g., "PX965180.1 Coxsackievirus...")
            acc_full = record.id.split()[0]
            acc_base = acc_full.split('.')[0]
            seqs.append((record, acc_full, acc_base))
    except Exception as e:
        print(f"Error reading FASTA: {e}")
        sys.exit(1)
    return seqs

def sort_sequences_and_metadata(
    fasta_file: str,
    tsv_file: str,
    field: str,
    output_fasta: str,
    output_tsv: str
):
    """Sort sequences and metadata by field."""
    print(f"Loading metadata from {tsv_file}...")
    meta = load_metadata(tsv_file)
    print(f"Loaded {len(meta):,} metadata records")

    print(f"Loading sequences from {fasta_file}...")
    seqs = load_sequences(fasta_file)
    print(f"Loaded {len(seqs):,} sequences")

    # Determine field type
    field_type = 'string'
    if field in ('taxid', 'length'):
        field_type = 'numeric'
    elif field == 'collection_date':
        field_type = 'date'
        # We'll sort by year extracted
        field = 'year_extracted'

    # Build sort keys
    def get_sort_key(item):
        record, acc_full, acc_base = item
        # Try to get metadata
        row = meta.get(acc_full) or meta.get(acc_base)
        if row is None:
            # No metadata: put at end
            return (1, '')
        if field == 'year_extracted':
            date_str = row.get('collection_date', '')
            year = parse_year(date_str)
            return (0, year)
        elif field in ('species', 'genus', 'organism', 'taxid', 'country', 'genotype', 'host'):
            val = row.get(field, '')
            return (0, val.lower() if isinstance(val, str) else val)
        elif field_type == 'numeric':
            try:
                val = row.get(field, '0')
                return (0, int(val))
            except (ValueError, TypeError):
                return (0, 0)
        else:
            val = row.get(field, '')
            if isinstance(val, str):
                val = val.lower()
            return (0, val)

    print(f"Sorting by field '{field}'...")
    sorted_seqs = sorted(seqs, key=get_sort_key)

    # Write sorted FASTA
    print(f"Writing sorted FASTA to {output_fasta}...")
    with open(output_fasta, 'w') as f:
        for record, _, _ in sorted_seqs:
            SeqIO.write(record, f, "fasta")

    # Write sorted metadata (in same order as sorted sequences)
    print(f"Writing sorted metadata to {output_tsv}...")
    # Get fieldnames from original TSV
    fieldnames = []
    try:
        with open(tsv_file, 'r', newline='') as f:
            reader = csv.DictReader(f, delimiter='\t')
            fieldnames = reader.fieldnames
    except Exception:
        fieldnames = ['accession', 'version', 'length', 'taxid', 'organism',
                     'species', 'genus', 'family', 'collection_date', 'country',
                     'isolate', 'strain', 'host', 'genotype', 'submitter', 'seq_tech']
    # Write header + rows in sorted order
    written = set()
    with open(output_tsv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for _, acc_full, acc_base in sorted_seqs:
            row = meta.get(acc_full) or meta.get(acc_base)
            if row:
                key = row.get('accession', '') or row.get('version', '')
                if key not in written:
                    writer.writerow(row)
                    written.add(key)

    print(f"\nDone. Sorted {len(sorted_seqs):,} sequences by '{field}'.")

def main():
    parser = argparse.ArgumentParser(
        description='Sort FASTA sequences and metadata TSV by a specified field.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Sort by taxid
  python sort_seqs.py -f seqs.fasta -m meta.tsv --field taxid -of sorted.fasta -om sorted.tsv

  # Sort by species
  python sort_seqs.py -f seqs.fasta -m meta.tsv --field species -of sorted.fasta -om sorted.tsv

  # Sort by year (extracted from collection_date)
  python sort_seqs.py -f seqs.fasta -m meta.tsv --field collection_date -of sorted.fasta -om sorted.tsv
        """
    )
    parser.add_argument('-f', '--fasta', required=True, help='Input FASTA file')
    parser.add_argument('-m', '--metadata', required=True, help='Input metadata TSV file')
    parser.add_argument('--field', required=True,
                        help='Field to sort by (taxid, organism, species, genus, collection_date, country, genotype, host, length, etc.)')
    parser.add_argument('-of', '--output-fasta', required=True, help='Output sorted FASTA file')
    parser.add_argument('-om', '--output-metadata', required=True, help='Output sorted metadata TSV file')
    args = parser.parse_args()

    if not os.path.exists(args.fasta):
        print(f"Error: FASTA file {args.fasta} not found")
        sys.exit(1)
    if not os.path.exists(args.metadata):
        print(f"Error: Metadata file {args.metadata} not found")
        sys.exit(1)

    sort_sequences_and_metadata(
        args.fasta,
        args.metadata,
        args.field,
        args.output_fasta,
        args.output_metadata
    )

if __name__ == '__main__':
    main()
