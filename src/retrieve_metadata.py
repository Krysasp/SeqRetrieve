#!/usr/bin/env python3
"""
Retrieve metadata for GenBank accession numbers.
Outputs TSV with: accession, version, taxid, species, genus, family,
collection_date, country, isolate, strain, submitter, genotype, length.
"""

import sys
import os
import argparse
import csv
import time
import random
from typing import List, Optional, Dict, Set
from Bio import Entrez, SeqIO
import signal

class MetadataFetcher:
    def __init__(self, email: str, max_retries: int = 5,
                 base_delay: float = 2.0, max_delay: float = 60.0):
        self.email = email
        self.max_retries = max_retries
        self.base_delay = base_delay
        self.max_delay = max_delay
        self.results = []
        self.failed = []
        Entrez.email = email
        Entrez.tool = "RetrieveSeqMetadata"
        signal.signal(signal.SIGINT, self.signal_handler)
        self.interrupted = False

    def signal_handler(self, sig, frame):
        print("\nInterrupted. Saving progress...")
        self.interrupted = True

    def exponential_backoff(self, attempt: int) -> float:
        delay = min(self.max_delay, self.base_delay * (2 ** attempt))
        jitter = random.uniform(0, delay * 0.1)
        return delay + jitter

    def read_accessions(self, input_file: str, column: str = None) -> List[str]:
        accessions = []
        try:
            with open(input_file, 'r') as f:
                sample = f.read(1024)
                f.seek(0)
                has_header = csv.Sniffer().has_header(sample)
                if has_header:
                    reader = csv.DictReader(f)
                    if column:
                        if column not in reader.fieldnames:
                            print(f"Error: Column '{column}' not found. Available: {reader.fieldnames}")
                            sys.exit(1)
                        accessions = [row[column].strip() for row in reader if row[column].strip()]
                    else:
                        first_col = reader.fieldnames[0]
                        print(f"Using column '{first_col}'")
                        accessions = [row[first_col].strip() for row in reader if row[first_col].strip()]
                else:
                    reader = csv.reader(f)
                    accessions = [row[0].strip() for row in reader if row and row[0].strip()]
        except Exception as e:
            print(f"Error reading input: {e}")
            sys.exit(1)
        # Deduplicate
        seen = set()
        unique = []
        for acc in accessions:
            if acc not in seen:
                seen.add(acc)
                unique.append(acc)
        print(f"Read {len(unique):,} unique accessions")
        return unique

    def extract_accessions_from_fasta(self, fasta_file: str) -> Set[str]:
        accs = set()
        try:
            for record in SeqIO.parse(fasta_file, "fasta"):
                acc = record.id.split()[0]
                accs.add(acc)
                accs.add(acc.split('.')[0])
        except Exception as e:
            print(f"Warning: Could not read FASTA: {e}")
        return accs

    def fetch_metadata_with_retry(self, accession: str, retry: int = 0):
        if retry >= self.max_retries:
            return None
        try:
            time.sleep(random.uniform(0.1, 0.3))
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
            record = next(SeqIO.parse(handle, "genbank"))
            handle.close()
            return self.parse_genbank_metadata(record, accession)
        except Exception as e:
            error_msg = str(e)
            if "429" in error_msg or "Too Many Requests" in error_msg:
                wait = self.exponential_backoff(retry)
                print(f"  Rate limited on {accession}, waiting {wait:.1f}s (attempt {retry+1})")
                time.sleep(wait)
                return self.fetch_metadata_with_retry(accession, retry + 1)
            elif "404" in error_msg or "Not Found" in error_msg:
                print(f"  Accession {accession} not found")
                return None
            else:
                if retry < self.max_retries - 1:
                    wait = self.exponential_backoff(retry) / 2
                    print(f"  Error on {accession}: {error_msg[:60]}, retrying in {wait:.1f}s")
                    time.sleep(wait)
                    return self.fetch_metadata_with_retry(accession, retry + 1)
                else:
                    print(f"  Failed {accession}: {error_msg[:100]}")
                    return None

    def parse_genbank_metadata(self, record, accession: str) -> Dict:
        data = {
            'accession': record.name,
            'version': record.id,
            'length': len(record.seq),
            'taxid': '',
            'species': '',
            'genus': '',
            'family': '',
            'collection_date': '',
            'country': '',
            'isolate': '',
            'strain': '',
            'submitter': '',
            'genotype': '',
            'organism': '',
            'host': '',
            'seq_tech': '',
        }
        # Organism
        if record.annotations.get('organism'):
            data['organism'] = record.annotations['organism']
            parts = record.annotations['organism'].split()
            if len(parts) >= 1:
                data['genus'] = parts[0]
            if len(parts) >= 2:
                data['species'] = ' '.join(parts[:2])
        # TaxID
        for xref in record.dbxrefs:
            if xref.startswith('taxon:'):
                data['taxid'] = xref.split(':')[1]
                break
        # Source feature qualifiers
        for feature in record.features:
            if feature.type == 'source':
                qual = feature.qualifiers
                data['collection_date'] = qual.get('collection_date', [''])[0]
                data['country'] = qual.get('country', [''])[0]
                data['isolate'] = qual.get('isolate', [''])[0]
                data['strain'] = qual.get('strain', [''])[0]
                data['host'] = qual.get('host', [''])[0]
                data['genotype'] = qual.get('genotype', [''])[0]
                if not data['genotype']:
                    notes = qual.get('note', [])
                    for n in notes:
                        if 'genotype' in n.lower():
                            data['genotype'] = n
                            break
                break
        # Submitter
        if 'references' in record.annotations:
            for ref in record.annotations['references']:
                if ref.authors:
                    data['submitter'] = ref.authors
                    break
        return data

    def fetch_batch_metadata(self, batch: List[str], batch_num: int, total_batches: int):
        results = []
        for attempt in range(self.max_retries):
            if self.interrupted:
                return results
            try:
                post_handle = Entrez.epost(db="nucleotide", id=",".join(batch))
                post_result = Entrez.read(post_handle)
                post_handle.close()
                time.sleep(0.5)
                fetch_handle = Entrez.efetch(
                    db="nucleotide",
                    rettype="gb",
                    retmode="text",
                    webenv=post_result["WebEnv"],
                    query_key=post_result["QueryKey"],
                    retmax=len(batch)
                )
                records = list(SeqIO.parse(fetch_handle, "genbank"))
                fetch_handle.close()
                for record in records:
                    data = self.parse_genbank_metadata(record, '')
                    results.append(data)
                # Check for missing
                got_accs = {d['accession'] for d in results}
                got_accs.update({d['version'] for d in results})
                failed = [acc for acc in batch if acc not in got_accs and acc.split('.')[0] not in got_accs]
                if not failed:
                    return results, failed
                # Retry individually for failed
                if attempt < self.max_retries - 1:
                    print(f"  Batch {batch_num}: {len(failed)} failed, retrying individually...")
                    for acc in failed:
                        d = self.fetch_metadata_with_retry(acc)
                        if d:
                            results.append(d)
                    got_accs = {d['accession'] for d in results}
                    got_accs.update({d['version'] for d in results})
                    failed = [acc for acc in batch if acc not in got_accs and acc.split('.')[0] not in got_accs]
                    if not failed:
                        return results, failed
                else:
                    self.failed.extend(failed)
                    return results, failed
            except Exception as e:
                error_msg = str(e)
                if "429" in error_msg or "Too Many Requests" in error_msg:
                    wait = self.exponential_backoff(attempt)
                    print(f"  Batch {batch_num} rate limited, waiting {wait:.1f}s")
                    time.sleep(wait)
                else:
                    print(f"  Batch {batch_num} error: {error_msg[:100]}")
                    if attempt < self.max_retries - 1:
                        time.sleep(self.exponential_backoff(attempt) / 2)
                    else:
                        self.failed.extend(batch)
                        return results, batch
        return results, []

    def fetch_all_metadata(self, accessions: List[str], batch_size: int = 50,
                           resume: bool = False, existing_tsv: str = None):
        if resume and existing_tsv and os.path.exists(existing_tsv):
            existing = self.read_existing_tsv(existing_tsv)
            to_fetch = [acc for acc in accessions if acc not in existing]
            print(f"Resume: {len(accessions) - len(to_fetch):,} already have metadata, fetching {len(to_fetch):,}")
            accessions = to_fetch
        if not accessions:
            print("Nothing to fetch.")
            return True
        batches = [accessions[i:i+batch_size] for i in range(0, len(accessions), batch_size)]
        total = len(batches)
        print(f"Fetching metadata: {len(accessions):,} accessions, {total} batches")
        start = time.time()
        for num, batch in enumerate(batches, 1):
            if self.interrupted:
                break
            elapsed = time.time() - start
            avg = elapsed / (num - 1) if num > 1 else 0
            remaining = total - num
            eta = avg * remaining
            print(f"\nBatch {num}/{total} | ETA: {eta/60:.1f} min | Progress: {len(self.results):,}/{len(accessions):,}")
            batch_results, _ = self.fetch_batch_metadata(batch, num, total)
            self.results.extend(batch_results)
            print(f"  Got {len(batch_results)} metadata records")
            if num < total:
                time.sleep(random.uniform(1.0, 2.0))
            if num % 10 == 0:
                self.save_intermediate_tsv(f"intermediate_meta_batch_{num}.tsv")
        return len(self.results) > 0

    def read_existing_tsv(self, tsv_file: str) -> Set[str]:
        accs = set()
        try:
            with open(tsv_file, 'r') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    accs.add(row.get('accession', ''))
                    accs.add(row.get('version', ''))
        except Exception as e:
            print(f"Warning: Could not read existing TSV: {e}")
        return accs

    def save_metadata_tsv(self, output_file: str, resume: bool = False):
        if not self.results:
            print("No metadata to save")
            return False
        fieldnames = ['accession', 'version', 'length', 'taxid', 'organism',
                      'species', 'genus', 'family', 'collection_date', 'country',
                      'isolate', 'strain', 'host', 'genotype', 'submitter', 'seq_tech']
        mode = 'a' if resume and os.path.exists(output_file) else 'w'
        try:
            with open(output_file, mode, newline='') as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
                if mode == 'w':
                    writer.writeheader()
                for d in self.results:
                    writer.writerow(d)
            print(f"\nSaved {len(self.results):,} metadata records to {output_file}")
            return True
        except Exception as e:
            print(f"Error saving: {e}")
            return False

    def save_intermediate_tsv(self, filename: str):
        if self.results:
            self.save_metadata_tsv(filename)

    def save_failed(self, filename: str = "failed_metadata.txt"):
        if self.failed:
            with open(filename, 'w') as f:
                for acc in set(self.failed):
                    f.write(acc + '\n')
            print(f"Failed accessions saved to {filename}")


def main():
    parser = argparse.ArgumentParser(
        description='Retrieve metadata for GenBank accession numbers',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  # Retrieve metadata for accessions in CSV
  python retrieve_metadata.py -i accessions.csv -o metadata.tsv -e your@email.com

  # Resume from existing TSV
  python retrieve_metadata.py -i accessions.csv -o metadata.tsv -e your@email.com -r metadata.tsv

  # Get metadata from existing FASTA file
  python retrieve_metadata.py -f sequences.fasta -o metadata.tsv -e your@email.com
        """
    )
    parser.add_argument('-i', '--input', help='Input CSV file with accession numbers')
    parser.add_argument('-f', '--fasta', help='Input FASTA file to extract accessions from')
    parser.add_argument('-o', '--output', required=True, help='Output TSV file')
    parser.add_argument('-e', '--email', required=True, help='Email for NCBI')
    parser.add_argument('-c', '--column', default=None, help='Column name for accessions')
    parser.add_argument('-r', '--resume', default=None, help='Resume from existing TSV')
    parser.add_argument('--batch-size', type=int, default=50, help='Batch size (default 50)')
    parser.add_argument('--max-retries', type=int, default=5, help='Max retries (default 5)')
    args = parser.parse_args()

    fetcher = MetadataFetcher(email=args.email, max_retries=args.max_retries)

    accessions = []
    if args.fasta:
        accs = fetcher.extract_accessions_from_fasta(args.fasta)
        accessions = list(accs)
        print(f"Extracted {len(accessions):,} accessions from FASTA")
    elif args.input:
        accessions = fetcher.read_accessions(args.input, args.column)
    else:
        print("Error: Need -i/--input or -f/--fasta")
        sys.exit(1)

    if not accessions:
        print("No accessions found")
        sys.exit(1)

    resume_mode = False
    if args.resume:
        if not os.path.exists(args.resume):
            print(f"Error: Resume file {args.resume} not found")
            sys.exit(1)
        resume_mode = True

    print(f"\n{'='*70}")
    print(f"METADATA RETRIEVER")
    print(f"{'='*70}")
    print(f"Input: {args.input or args.fasta}")
    print(f"Output: {args.output}")
    print(f"Total accessions: {len(accessions):,}")
    if resume_mode:
        print(f"Resume file: {args.resume}")
    print(f"{'='*70}")

    success = fetcher.fetch_all_metadata(accessions, batch_size=args.batch_size,
                                        resume=resume_mode, existing_tsv=args.resume)
    if success:
        fetcher.save_metadata_tsv(args.output, resume=resume_mode)
        fetcher.save_failed()
        print(f"\nDone. Retrieved {len(fetcher.results):,} metadata records.")
        if fetcher.failed:
            print(f"Failed: {len(set(fetcher.failed)):,}")
            sys.exit(1)
    else:
        print("Failed to retrieve metadata")
        sys.exit(1)


if __name__ == '__main__':
    main()
