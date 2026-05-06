#!/usr/bin/env python3
"""
Auto-align sequences grouped by a specified field (taxid, organism, species, genus, year, etc.).
Uses external alignment tools: MUSCLE, Clustal Omega, or MAFFT.
"""

import sys
import os
import argparse
import csv
import subprocess
from Bio import SeqIO
from typing import Dict, List, Tuple
import re
import tempfile
import shutil

def parse_year(date_str: str) -> str:
    if not date_str:
        return 'unknown'
    match = re.search(r'(\d{4})', date_str)
    if match:
        return match.group(1)
    return 'unknown'

def load_metadata(tsv_file: str) -> Dict[str, dict]:
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
                    base = ver.split('.')[0]
                    if base not in meta:
                        meta[base] = row
    except Exception as e:
        print(f"Error reading metadata: {e}")
        sys.exit(1)
    return meta

def group_sequences(fasta_file: str, meta: Dict, field: str) -> Dict[str, List[Tuple[str, object]]]:
    """Group sequences by field. Returns dict: group_key -> list of (record_id, record)."""
    groups = {}
    unknown_key = 'unknown'
    # Determine if field needs special handling
    if field == 'collection_date':
        # Group by year
        field_key = 'year'
    else:
        field_key = field

    for record in SeqIO.parse(fasta_file, "fasta"):
        acc_full = record.id.split()[0]
        acc_base = acc_full.split('.')[0]
        row = meta.get(acc_full) or meta.get(acc_base)
        group = None
        if row:
            if field == 'collection_date':
                date_str = row.get('collection_date', '')
                group = parse_year(date_str)
            elif field in ('species', 'genus', 'organism', 'taxid', 'country', 'genotype', 'host'):
                group = row.get(field, '').strip()
            else:
                group = row.get(field, '').strip()
        if not group:
            group = unknown_key
        # Sanitize group name for filename
        safe_group = "".join(c if c.isalnum() or c in '_- ' else '_' for c in group).strip()
        if not safe_group:
            safe_group = unknown_key
        groups.setdefault(safe_group, []).append(record)
    return groups

def find_alignment_tool(preferred: str = None):
    """Find available alignment tool."""
    tools = []
    if preferred:
        tools.append(preferred)
    tools.extend(['muscle', 'clustalo', 'mafft'])
    for tool in tools:
        if shutil.which(tool):
            return tool
    return None

def align_group(group_name: str, records: List, tool: str, output_dir: str, temp_dir: str) -> str:
    """Align a group of sequences. Returns path to aligned file."""
    if len(records) < 2:
        # No alignment needed for single sequence
        out_file = os.path.join(output_dir, f"{group_name}.fasta")
        with open(out_file, 'w') as f:
            SeqIO.write(records, f, "fasta")
        print(f"  Group '{group_name}': single sequence, copied (no alignment)")
        return out_file

    # Write input fasta
    input_fasta = os.path.join(temp_dir, f"input_{group_name}.fasta")
    with open(input_fasta, 'w') as f:
        SeqIO.write(records, f, "fasta")

    out_file = os.path.join(output_dir, f"{group_name}.aligned.fasta")
    out_file_clustal = os.path.join(temp_dir, f"output_{group_name}.clustal")

    try:
        if tool == 'muscle':
            subprocess.run(['muscle', '-in', input_fasta, '-out', out_file],
                          check=True, capture_output=True)
        elif tool == 'clustalo':
            subprocess.run(['clustalo', '-i', input_fasta, '-o', out_file, '--outfmt=fasta'],
                          check=True, capture_output=True)
        elif tool == 'mafft':
            with open(out_file, 'w') as outf:
                subprocess.run(['mafft', '--auto', input_fasta], stdout=outf,
                              check=True, capture_output=False)
        else:
            print(f"  Unsupported tool: {tool}")
            return None
        print(f"  Group '{group_name}': aligned {len(records)} sequences")
        return out_file
    except subprocess.CalledProcessError as e:
        print(f"  Error aligning group '{group_name}': {e}")
        # Fallback: copy unaligned
        with open(out_file, 'w') as f:
            SeqIO.write(records, f, "fasta")
        return out_file

def main():
    parser = argparse.ArgumentParser(
        description='Auto-align sequences grouped by a specified field.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Align by species using MUSCLE
  python align_seqs.py -f seqs.fasta -m meta.tsv --field species -o aligned/

  # Align by taxid using Clustal Omega
  python align_seqs.py -f seqs.fasta -m meta.tsv --field taxid -o aligned/ --tool clustalo

  # Align by year (from collection_date)
  python align_seqs.py -f seqs.fasta -m meta.tsv --field collection_date -o aligned/
        """
    )
    parser.add_argument('-f', '--fasta', required=True, help='Input FASTA file')
    parser.add_argument('-m', '--metadata', required=True, help='Metadata TSV file')
    parser.add_argument('--field', required=True,
                        help='Field to group by (taxid, organism, species, genus, collection_date, country, etc.)')
    parser.add_argument('-o', '--output-dir', required=True, help='Output directory for aligned files')
    parser.add_argument('--tool', default=None,
                        help='Alignment tool (muscle, clustalo, mafft). Auto-detect if not specified.')
    parser.add_argument('--combined-output', default=None,
                        help='Optional: combined aligned FASTA file for all groups')
    args = parser.parse_args()

    if not os.path.exists(args.fasta):
        print(f"Error: FASTA file {args.fasta} not found")
        sys.exit(1)
    if not os.path.exists(args.metadata):
        print(f"Error: Metadata file {args.metadata} not found")
        sys.exit(1)

    # Find alignment tool
    tool = args.tool
    if not tool:
        tool = find_alignment_tool()
        if not tool:
            print("Error: No alignment tool found. Install muscle, clustalo, or mafft.")
            sys.exit(1)
    else:
        if not shutil.which(tool):
            print(f"Error: Tool '{tool}' not found in PATH")
            sys.exit(1)

    print(f"Using alignment tool: {tool}")

    # Create output dir
    os.makedirs(args.output_dir, exist_ok=True)

    # Load metadata
    print(f"Loading metadata from {args.metadata}...")
    meta = load_metadata(args.metadata)
    print(f"Loaded {len(meta):,} metadata records")

    # Group sequences
    print(f"Grouping sequences by '{args.field}'...")
    groups = group_sequences(args.fasta, meta, args.field)
    print(f"Found {len(groups)} groups:")
    for grp, recs in sorted(groups.items()):
        print(f"  {grp}: {len(recs)} sequences")

    # Align each group
    print(f"\nAligning {len(groups)} groups...")
    temp_dir = tempfile.mkdtemp(prefix="align_")
    aligned_files = []
    try:
        for grp, recs in sorted(groups.items()):
            out = align_group(grp, recs, tool, args.output_dir, temp_dir)
            if out:
                aligned_files.append(out)
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)

    # Combined output
    if args.combined_output:
        print(f"\nWriting combined aligned sequences to {args.combined_output}...")
        with open(args.combined_output, 'w') as outf:
            for af in aligned_files:
                for record in SeqIO.parse(af, "fasta"):
                    SeqIO.write(record, outf, "fasta")
        print(f"Combined file written with {sum(len(list(SeqIO.parse(af, 'fasta'))) for af in aligned_files)} sequences")

    print(f"\nDone. Aligned files in {args.output_dir}")

if __name__ == '__main__':
    main()
