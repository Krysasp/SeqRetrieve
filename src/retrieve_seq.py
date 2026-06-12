#!/usr/bin/env python3
"""
Optimized GenBank sequence retriever with rate limiting and intelligent retry strategies
Handles 20,000+ sequences without triggering HTTP 429 errors
Supports resume mode to continue from a previous run
"""

import sys
import os
import argparse
import csv
import time
import random
from typing import List, Optional, Set, Dict, Tuple
from Bio import Entrez, SeqIO
import signal
from datetime import datetime

class SmartGenBankFetcher:
    def __init__(self, email: str, max_retries: int = 5, 
                 base_delay: float = 2.0, max_delay: float = 60.0):
        """
        Initialize the smart fetcher with exponential backoff.
        
        Args:
            email: Email address for NCBI
            max_retries: Maximum number of retry attempts
            base_delay: Initial delay between requests (seconds)
            max_delay: Maximum delay for exponential backoff
        """
        self.email = email
        self.max_retries = max_retries
        self.base_delay = base_delay
        self.max_delay = max_delay
        self.results = []
        self.failed_accessions = []
        self.successful_accessions = set()
        
        # Configure Entrez
        Entrez.email = email
        Entrez.tool = "OptimizedGenBankRetriever"
        
        # Set up signal handler for graceful interruption
        signal.signal(signal.SIGINT, self.signal_handler)
        self.interrupted = False
        
    def signal_handler(self, sig, frame):
        """Handle Ctrl+C gracefully."""
        print("\n\nInterrupted by user. Saving current progress...")
        self.interrupted = True
    
    def read_accessions_from_csv(self, csv_file: str, column_name: str = None) -> List[str]:
        """Read accession numbers from CSV file."""
        accessions = []

        try:
            with open(csv_file, 'r') as f:
                # Check for header
                sample = f.read(1024)
                f.seek(0)
                has_header = csv.Sniffer().has_header(sample)

                if has_header:
                    reader = csv.DictReader(f)
                    if column_name:
                        if column_name not in reader.fieldnames:
                            print(f"Error: Column '{column_name}' not found")
                            print(f"Available columns: {', '.join(reader.fieldnames)}")
                            sys.exit(1)
                        accessions = [row[column_name].strip() for row in reader if row[column_name].strip()]
                    else:
                        first_col = reader.fieldnames[0]
                        print(f"Using column '{first_col}' for accession numbers")
                        accessions = [row[first_col].strip() for row in reader if row[first_col].strip()]
                else:
                    reader = csv.reader(f)
                    accessions = [row[0].strip() for row in reader if row and row[0].strip()]

        except Exception as e:
            print(f"Error reading CSV: {e}")
            sys.exit(1)

        # Remove duplicates while preserving order
        seen = set()
        unique_accessions = []
        for acc in accessions:
            if acc not in seen:
                seen.add(acc)
                unique_accessions.append(acc)

        print(f"Read {len(unique_accessions):,} unique accession numbers from {csv_file}")
        return unique_accessions

    def extract_accessions_from_fasta(self, fasta_file: str) -> Set[str]:
        """Extract accession numbers (with and without version) from an existing FASTA file."""
        accessions = set()
        try:
            for record in SeqIO.parse(fasta_file, "fasta"):
                # Extract accession from record.id (e.g., "PX965180.1" -> "PX965180", "PX965180.1")
                acc_with_version = record.id.split()[0]
                accessions.add(acc_with_version)
                # Also add without version
                acc_no_version = acc_with_version.split('.')[0]
                accessions.add(acc_no_version)
            print(f"Found {len(accessions):,} accession variants in existing FASTA file")
        except Exception as e:
            print(f"Warning: Could not read existing FASTA file: {e}")
        return accessions
    
    def compile_accessions_from_output_dir(self, output_dir: str) -> Set[str]:
        """
        Compile accession IDs from all FASTA files in the output directory.
        Prompts user to confirm compilation.
        
        Args:
            output_dir: Directory containing FASTA files
            
        Returns:
            Set of accession numbers from all FASTA files
        """
        if not os.path.isdir(output_dir):
            print(f"Warning: Output directory '{output_dir}' not found")
            return set()
        
        # Find all FASTA files
        fasta_files = []
        for f in os.listdir(output_dir):
            if f.endswith('.fasta') or f.endswith('.fa'):
                fasta_files.append(os.path.join(output_dir, f))
        
        if not fasta_files:
            print(f"No FASTA files found in {output_dir}")
            return set()
        
        print(f"\nFound {len(fasta_files)} FASTA file(s) in {output_dir}:")
        for f in fasta_files:
            file_size = os.path.getsize(f) / (1024 * 1024)
            print(f"  - {os.path.basename(f)} ({file_size:.1f} MB)")
        
        # Prompt user
        response = input("\nCompile accession IDs from these files? (y/n): ")
        if response.lower() != 'y':
            print("Skipping compilation")
            return set()
        
        # Compile accessions from all files
        all_accessions = set()
        total_sequences = 0
        
        for fasta_file in fasta_files:
            try:
                for record in SeqIO.parse(fasta_file, "fasta"):
                    acc_with_version = record.id.split()[0]
                    all_accessions.add(acc_with_version)
                    acc_no_version = acc_with_version.split('.')[0]
                    all_accessions.add(acc_no_version)
                    total_sequences += 1
                print(f"  Loaded {len(all_accessions):,} accessions from {os.path.basename(fasta_file)}")
            except Exception as e:
                print(f"  Warning: Could not read {os.path.basename(fasta_file)}: {e}")
        
        print(f"\nTotal compiled accessions: {len(all_accessions):,} (from {total_sequences} sequences)")
        return all_accessions
    
    def load_existing_sequences(self, output_dir: str, target_accessions: Set[str]) -> List[object]:
        """
        Load existing sequences from FASTA files in the output directory.
        Deduplicates sequences by base accession ID (keeps first occurrence).
        
        Args:
            output_dir: Directory containing FASTA files
            target_accessions: Set of accession IDs to load (base accessions without version)
            
        Returns:
            List of SeqRecord objects, deduplicated by base accession ID
        """
        existing_records = []
        seen_accessions = set()
        
        # Find all FASTA files
        fasta_files = []
        for f in os.listdir(output_dir):
            if f.endswith('.fasta') or f.endswith('.fa'):
                # Exclude intermediate checkpoint files and current run files
                if '_run_' not in f or not f.startswith(datetime.now().strftime('%Y-%m-%d')):
                    fasta_files.append(os.path.join(output_dir, f))
        
        for fasta_file in sorted(fasta_files):
            try:
                for record in SeqIO.parse(fasta_file, "fasta"):
                    acc_with_version = record.id.split()[0]
                    acc_no_version = acc_with_version.split('.')[0]
                    
                    # Only include if in target accessions and not already seen
                    if acc_no_version in target_accessions and acc_no_version not in seen_accessions:
                        existing_records.append(record)
                        seen_accessions.add(acc_no_version)
            except Exception as e:
                print(f"  Warning: Could not read {os.path.basename(fasta_file)}: {e}")
        
        print(f"  Loaded {len(existing_records):,} unique sequences (deduplicated by accession ID)")
        return existing_records
    
    def exponential_backoff(self, attempt: int) -> float:
        """Calculate delay with exponential backoff and jitter."""
        delay = min(self.max_delay, self.base_delay * (2 ** attempt))
        # Add random jitter to avoid thundering herd
        jitter = random.uniform(0, delay * 0.1)
        return delay + jitter
    
    def fetch_with_retry(self, accession: str, retry_count: int = 0) -> Optional[object]:
        """
        Fetch a single sequence with intelligent retry logic.
        """
        if retry_count >= self.max_retries:
            return None
            
        try:
            # Add small random delay before request to spread load
            time.sleep(random.uniform(0.1, 0.3))
            
            handle = Entrez.efetch(
                db="nucleotide",
                id=accession,
                rettype="fasta",
                retmode="text"
            )
            record = SeqIO.read(handle, "fasta")
            handle.close()
            
            # Success - reset retry count for this accession
            return record
            
        except Exception as e:
            error_msg = str(e)
            
            # Handle rate limiting specially
            if "429" in error_msg or "Too Many Requests" in error_msg:
                wait_time = self.exponential_backoff(retry_count)
                print(f"  Rate limited on {accession}, waiting {wait_time:.1f}s (attempt {retry_count + 1}/{self.max_retries})")
                time.sleep(wait_time)
                return self.fetch_with_retry(accession, retry_count + 1)
                
            elif "502" in error_msg or "Bad Gateway" in error_msg:
                wait_time = self.exponential_backoff(retry_count) / 2
                print(f"  Server error on {accession}, waiting {wait_time:.1f}s (attempt {retry_count + 1}/{self.max_retries})")
                time.sleep(wait_time)
                return self.fetch_with_retry(accession, retry_count + 1)
                
            elif "404" in error_msg or "Not Found" in error_msg:
                # Accession doesn't exist, don't retry
                print(f"  Accession {accession} not found")
                return None
                
            else:
                # Other errors, retry with backoff
                if retry_count < self.max_retries - 1:
                    wait_time = self.exponential_backoff(retry_count) / 2
                    print(f"  Error on {accession}: {error_msg[:50]}, retrying in {wait_time:.1f}s")
                    time.sleep(wait_time)
                    return self.fetch_with_retry(accession, retry_count + 1)
                else:
                    print(f"  Failed to fetch {accession}: {error_msg[:100]}")
                    return None
    
    def fetch_batch_smart(self, batch: List[str], batch_num: int, total_batches: int) -> tuple:
        """
        Fetch a batch of sequences using NCBI's batch processing with rate limiting.
        Waits until all sequences are retrieved or max retries exceeded for the batch.
        """
        batch_results = []
        batch_failed = list(batch)  # Start with all, remove as they succeed

        for attempt in range(self.max_retries):
            if self.interrupted:
                return batch_results, batch_failed

            if not batch_failed:
                return batch_results, []

            # Longer initial delay for batch requests
            time.sleep(random.uniform(1.0, 2.0))

            try:
                # Use epost for batch retrieval
                post_handle = Entrez.epost(db="nucleotide", id=",".join(batch_failed))
                post_result = Entrez.read(post_handle)
                post_handle.close()

                # Small delay between post and fetch
                time.sleep(0.5)

                fetch_handle = Entrez.efetch(
                    db="nucleotide",
                    rettype="fasta",
                    retmode="text",
                    webenv=post_result["WebEnv"],
                    query_key=post_result["QueryKey"],
                    retmax=len(batch_failed)
                )

                # Parse results
                records = list(SeqIO.parse(fetch_handle, "fasta"))
                fetch_handle.close()

                # Map successful accessions
                successful_in_batch = set()
                for record in records:
                    batch_results.append(record)
                    # Extract accession from record ID
                    for acc in batch_failed:
                        if acc in record.id or acc.split('.')[0] in record.id:
                            successful_in_batch.add(acc)
                            break

                # Update failed list
                batch_failed = [acc for acc in batch_failed if acc not in successful_in_batch]

                if not batch_failed:
                    return batch_results, []

                # Some still failed, retry individually
                if attempt < self.max_retries - 1:
                    print(f"  Batch {batch_num}: {len(batch_failed)} failed, retrying individually...")
                    still_failed = []
                    for acc in batch_failed:
                        record = self.fetch_with_retry(acc)
                        if record:
                            batch_results.append(record)
                        else:
                            still_failed.append(acc)
                    batch_failed = still_failed

                    if not batch_failed:
                        return batch_results, []

            except Exception as e:
                error_msg = str(e)

                if "429" in error_msg or "Too Many Requests" in error_msg:
                    wait_time = self.exponential_backoff(attempt)
                    print(f"  Batch {batch_num} rate limited, waiting {wait_time:.1f}s (attempt {attempt + 1}/{self.max_retries})")
                    time.sleep(wait_time)

                elif "502" in error_msg or "Bad Gateway" in error_msg:
                    wait_time = self.exponential_backoff(attempt) / 2
                    print(f"  Batch {batch_num} server error, waiting {wait_time:.1f}s (attempt {attempt + 1}/{self.max_retries})")
                    time.sleep(wait_time)

                else:
                    print(f"  Batch {batch_num} error: {error_msg[:100]}")
                    if attempt < self.max_retries - 1:
                        wait_time = self.exponential_backoff(attempt) / 2
                        time.sleep(wait_time)

        return batch_results, batch_failed
    
    def fetch_batch_individual_fallback(self, batch: List[str], batch_num: int) -> tuple:
        """
        Fallback method: fetch each accession individually with rate limiting.
        """
        batch_results = []
        batch_failed = []
        
        print(f"  Batches {batch_num} - Processing {len(batch)} accessions individually...")
        
        for i, acc in enumerate(batch, 1):
            if self.interrupted:
                break
                
            # Show progress for individual fetches
            if i % 50 == 0:
                print(f"    Individual fetch progress: {i}/{len(batch)}", end='\r')
            
            record = self.fetch_with_retry(acc)
            if record:
                batch_results.append(record)
            else:
                batch_failed.append(acc)
            
            # Ensure we don't exceed rate limits
            time.sleep(random.uniform(0.2, 0.4))
        
        if len(batch) > 50:
            print()  # New line after progress indicator
            
        return batch_results, batch_failed
    
    def fetch_sequences(self, accessions: List[str], batch_size: int = 50, resume: bool = False, existing_fasta: str = None, date_str: str = None) -> bool:
        """
        Main fetch method with adaptive batching and rate limiting.

        Args:
            accessions: List of accession numbers
            batch_size: Size of batches (smaller is safer for rate limits)
            resume: If True, filter out already-retrieved accessions
            existing_fasta: Path to existing FASTA file for resume mode
            date_str: Date string for naming intermediate files (YYYY-MM-DD format)
        """
        # Handle resume mode
        if resume and existing_fasta:
            existing_accessions = self.extract_accessions_from_fasta(existing_fasta)
            accessions_to_fetch = [acc for acc in accessions if acc not in existing_accessions]
            skipped = len(accessions) - len(accessions_to_fetch)
            print(f"\nResume mode: {skipped:,} sequences already retrieved, {len(accessions_to_fetch):,} to fetch")
            accessions = accessions_to_fetch

        if not accessions:
            print("All sequences already retrieved. Nothing to do.")
            return True

        # Create batches
        batches = [accessions[i:i + batch_size] for i in range(0, len(accessions), batch_size)]
        total_batches = len(batches)

        print(f"\nStarting optimized retrieval with rate limiting strategy...")
        print(f"Total sequences to fetch: {len(accessions):,}")
        print(f"Batch size: {batch_size} (smaller batches to avoid rate limiting)")
        print(f"Total batches: {total_batches:,}")
        print("-" * 70)

        start_time = time.time()

        # Process batches sequentially with delays
        for batch_num, batch in enumerate(batches, 1):
            if self.interrupted:
                print("\nStopping due to user interruption...")
                break

            # Calculate ETA
            elapsed = time.time() - start_time
            avg_time_per_batch = elapsed / (batch_num - 1) if batch_num > 1 else 0
            remaining_batches = total_batches - batch_num
            eta_seconds = avg_time_per_batch * remaining_batches

            print(f"\nBatch {batch_num}/{total_batches} | "
                  f"ETA: {eta_seconds/60:.1f} min | "
                  f"Progress: {len(self.results):,}/{len(accessions):,} sequences")

            # Fetch batch (waits until all retrieved or max retries)
            batch_results, batch_failed = self.fetch_batch_smart(batch, batch_num, total_batches)

            # Update results
            self.results.extend(batch_results)
            self.failed_accessions.extend(batch_failed)

            # Track successful accessions
            for record in batch_results:
                for acc in batch:
                    if acc in record.id or acc.split('.')[0] in record.id:
                        self.successful_accessions.add(acc)
                        break

            # Show batch summary
            success_count = len(batch_results)
            fail_count = len(batch_failed)
            cumulative_success_rate = (len(self.results) / len(accessions)) * 100 if accessions else 0

            print(f"  Batch result: ✓ {success_count} | ✗ {fail_count} | "
                  f"Cumulative success: {len(self.results):,}/{len(accessions):,} ({cumulative_success_rate:.1f}%)")

            # Adaptive delay between batches based on failure rate
            if fail_count > len(batch) * 0.5:  # High failure rate
                delay = random.uniform(5.0, 10.0)
                print(f"  High failure rate detected, waiting {delay:.1f}s before next batch...")
                time.sleep(delay)
            else:
                # Normal delay
                delay = random.uniform(1.0, 2.0)
                time.sleep(delay)

            # Save intermediate results every 10 batches
            if batch_num % 10 == 0:
                run_num = batch_num // 10
                # Use provided date_str or generate one
                checkpoint_date = date_str if date_str else datetime.now().strftime('%Y-%m-%d')
                self.save_intermediate_results(f"{checkpoint_date}_run_{run_num}.fasta", checkpoint_date)

        # Calculate final statistics
        elapsed = time.time() - start_time
        success_rate = (len(self.results) / len(accessions)) * 100 if accessions else 0

        print(f"\n{'='*70}")
        print(f"RETRIEVAL COMPLETE")
        print(f"{'='*70}")
        print(f"Time elapsed: {elapsed:.1f} seconds ({elapsed/60:.1f} minutes)")
        print(f"Success rate: {success_rate:.1f}% ({len(self.results):,}/{len(accessions):,})")

        return len(self.results) > 0
    
    def save_results(self, output_file: str, resume: bool = False):
        """Save all successful sequences to a FASTA file. Append if resume mode."""
        if not self.results:
            print("No sequences to save")
            return False

        try:
            mode = 'a' if resume and os.path.exists(output_file) else 'w'
            with open(output_file, mode) as f:
                SeqIO.write(self.results, f, "fasta")

            file_size = os.path.getsize(output_file) / (1024 * 1024)  # MB
            action = "Appended" if mode == 'a' else "Saved"
            print(f"\n✓ {action} {len(self.results):,} sequences to {output_file}")
            print(f"  File size: {file_size:.1f} MB")
            return True

        except Exception as e:
            print(f"Error saving to {output_file}: {e}")
            return False
    
    def save_intermediate_results(self, filename: str, date_str: str = None):
        """Save intermediate results as checkpoint with date-based naming."""
        if self.results:
            try:
                if date_str:
                    # Use date-based naming for intermediate files
                    run_num = int(filename.split('_run_')[1].split('.')[0])
                    full_filename = f"{filename}"
                else:
                    full_filename = filename
                    
                with open(full_filename, 'w') as f:
                    SeqIO.write(self.results, f, "fasta")
                print(f"  ✓ Intermediate checkpoint saved to {full_filename}")
            except Exception as e:
                print(f"  ⚠ Could not save checkpoint: {e}")
    
    def save_failed_accessions(self, filename: str = "failed_accessions.txt"):
        """Save list of failed accessions for retry."""
        if self.failed_accessions:
            # Remove duplicates
            unique_failed = list(set(self.failed_accessions))
            with open(filename, 'w') as f:
                f.write('\n'.join(unique_failed))
            print(f"\nFailed accessions saved to {filename} ({len(unique_failed):,} sequences)")
            
            # Create a retry script
            retry_script = f"""#!/bin/bash
# Retry script for failed accessions
# Generated: {datetime.now()}
python3 {sys.argv[0]} -i {filename} -o retry_sequences.fasta -e {self.email} --batch-size 30
"""
            with open("retry_failed.sh", 'w') as f:
                f.write(retry_script)
            print(f"Retry script created: retry_failed.sh")
    
    def display_sequence_stats(self):
        """Display detailed statistics about retrieved sequences."""
        if not self.results:
            return
        
        print("\n" + "="*70)
        print("SEQUENCE STATISTICS")
        print("="*70)
        
        # Calculate statistics
        lengths = [len(seq.seq) for seq in self.results]
        total_bases = sum(lengths)
        min_len = min(lengths)
        max_len = max(lengths)
        avg_len = total_bases / len(lengths)
        
        # GC content (sample first 1000 sequences or less)
        sample_size = min(1000, len(self.results))
        sample_sequences = self.results[:sample_size]
        gc_contents = []
        for seq in sample_sequences:
            seq_str = str(seq.seq).upper()
            gc = (seq_str.count('G') + seq_str.count('C')) / len(seq_str) * 100
            gc_contents.append(gc)
        
        avg_gc = sum(gc_contents) / len(gc_contents) if gc_contents else 0
        
        print(f"Total sequences: {len(self.results):,}")
        print(f"Total bases: {total_bases:,} bp")
        print(f"Length range: {min_len:,} - {max_len:,} bp")
        print(f"Average length: {avg_len:.1f} bp")
        print(f"Average GC content: {avg_gc:.1f}% (based on {sample_size} sequences)")
        
        # Length distribution
        print("\nLength distribution:")
        ranges = [(0, 500), (500, 1000), (1000, 5000), (5000, 10000), (10000, float('inf'))]
        for low, high in ranges:
            if high == float('inf'):
                count = sum(1 for l in lengths if l >= low)
                label = f">{low:,}"
            else:
                count = sum(1 for l in lengths if low <= l < high)
                label = f"{low:,}-{high:,}"
            percent = (count / len(lengths)) * 100
            bar = '█' * int(percent / 2)
            print(f"  {label:>12}: {count:6,} ({percent:5.1f}%) {bar}")

def main():
    parser = argparse.ArgumentParser(
        description='Optimized GenBank sequence retriever with rate limiting and resume support',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
STRATEGIES TO AVOID RATE LIMITING:
  • Small batch sizes (30-50) to avoid triggering NCBI limits
  • Exponential backoff with jitter for retries
  • Adaptive delays between batches based on failure rates
  • Sequential processing (no parallel requests) to prevent 429 errors
  • Random delays to spread request timing
  • Waits for all sequences in a batch before proceeding to the next
  • Resume mode to continue from a previous run

RECOMMENDED USAGE:
  # Standard retrieval (safe for 20K+ sequences)
  python genbank_fetcher.py -i accessions.csv -o all_sequences.fasta -e your@email.com

  # More conservative (if still getting 429 errors)
  python genbank_fetcher.py -i accessions.csv -o sequences.fasta -e your@email.com --batch-size 30 --delay 5

  # Resume after interruption (continues from existing FASTA)
  python genbank_fetcher.py -i accessions.csv -o record_accession.fasta -e your@email.com -r record_accession.fasta

  # Resume with custom output
  python genbank_fetcher.py -i accessions.csv -o sequences.fasta -e your@email.com --resume sequences.fasta
        """
    )

    parser.add_argument('-i', '--input', required=True,
                       help='Input CSV file containing accession numbers')

    parser.add_argument('-o', '--output', required=True,
                       help='Output FASTA file for all sequences')

    parser.add_argument('-e', '--email', required=True,
                       help='Your email address (required by NCBI)')

    parser.add_argument('-c', '--column', default=None,
                       help='Column name containing accession numbers')

    parser.add_argument('-r', '--resume', default=None,
                       help='Resume mode: specify existing FASTA file to continue from')

    parser.add_argument('--batch-size', type=int, default=50,
                       help='Batch size for retrieval (default: 50, range: 20-100)')

    parser.add_argument('--base-delay', type=float, default=2.0,
                       help='Base delay between requests in seconds (default: 2.0)')

    parser.add_argument('--max-delay', type=float, default=60.0,
                       help='Maximum delay for exponential backoff (default: 60.0)')

    parser.add_argument('--max-retries', type=int, default=5,
                       help='Maximum retry attempts per batch (default: 5)')

    parser.add_argument('--test', action='store_true',
                       help='Test mode: only fetch first 100 sequences')

    parser.add_argument('--compile-existing', action='store_true',
                       help='Compile accession IDs from existing FASTA files in output directory before fetching')

    parser.add_argument('--output-dir', default=None,
                       help='Directory containing existing FASTA files for compilation (default: output/)')

    args = parser.parse_args()

    # Validate batch size
    if args.batch_size < 20:
        print("Warning: Batch size too small may be inefficient. Increasing to 20.")
        args.batch_size = 20
    elif args.batch_size > 100:
        print("Warning: Large batch sizes (>100) may trigger rate limiting. Consider using 30-50.")

    # Create fetcher
    fetcher = SmartGenBankFetcher(
        email=args.email,
        max_retries=args.max_retries,
        base_delay=args.base_delay,
        max_delay=args.max_delay
    )

    # Compile accessions from existing FASTA files if requested
    compiled_accessions = set()
    if args.compile_existing:
        output_dir = args.output_dir if args.output_dir else 'output'
        print(f"\n{'='*70}")
        print(f"COMPILING EXISTING FASTA FILES")
        print(f"{'='*70}")
        compiled_accessions = fetcher.compile_accessions_from_output_dir(output_dir)
        if not compiled_accessions:
            print("No accessions compiled, proceeding with CSV only")

    # Read accessions from CSV
    csv_accessions = fetcher.read_accessions_from_csv(args.input, args.column)

    if not csv_accessions:
        print("Error: No accession numbers found in CSV")
        sys.exit(1)

    # Combine accessions if compilation was performed
    if args.compile_existing and compiled_accessions:
        # Pre-filter: only fetch accessions NOT already in compiled set
        new_accessions = [acc for acc in csv_accessions if acc not in compiled_accessions]
        print(f"\n{'='*70}")
        print(f"ACCESSION SET ANALYSIS")
        print(f"{'='*70}")
        print(f"From existing FASTA (already downloaded): {len(compiled_accessions):,} accessions")
        print(f"From CSV (new to fetch): {len(new_accessions):,} accessions")
        print(f"Total unique in final output: {len(compiled_accessions) + len(new_accessions):,} accessions")
        
        # Store for later use: we'll load existing sequences and append new ones
        accessions_to_fetch = new_accessions
        accessions = csv_accessions  # Keep original for reference
    else:
        accessions_to_fetch = csv_accessions
        accessions = csv_accessions

    if not accessions:
        print("Error: No accession numbers to fetch")
        sys.exit(1)

    # Test mode
    if args.test:
        accessions = accessions[:100]
        print(f"\n*** TEST MODE: Limited to first 100 sequences ***")

    # Resume mode
    resume_mode = False
    if args.resume:
        if not os.path.exists(args.resume):
            print(f"Error: Resume file '{args.resume}' not found")
            sys.exit(1)
        resume_mode = True
        print(f"\n{'='*70}")
        print(f"RESUME MODE: Continuing from {args.resume}")
        print(f"{'='*70}")

    # Determine output filename and date string
    final_output = args.output
    date_str = datetime.now().strftime('%Y-%m-%d')
    
    if args.compile_existing and compiled_accessions:
        # Generate date-based filename with running number
        base_name = os.path.splitext(os.path.basename(args.output))[0]
        # Find existing run numbers in the final output directory to determine next run number
        # The final output will be created in the same directory as args.output
        output_dir = os.path.dirname(args.output)
        if not output_dir:
            output_dir = '.'  # Current directory
        existing_runs = []
        if os.path.exists(output_dir):
            for f in os.listdir(output_dir):
                if f.startswith(f"{date_str}_{base_name}") and f.endswith('.fasta'):
                    # Extract run number from filename (format: date_name_runN.fasta)
                    try:
                        suffix = f.replace(f"{date_str}_{base_name}_", "").replace('.fasta', "")
                        if suffix.startswith('run') and suffix[3:].isdigit():
                            existing_runs.append(int(suffix[3:]))
                    except:
                        pass
        run_num = max(existing_runs, default=0) + 1
        # Preserve the directory path from args.output
        if output_dir and output_dir != '.':
            final_output = f"{output_dir}/{date_str}_{base_name}_run{run_num}.fasta"
        else:
            final_output = f"{date_str}_{base_name}_run{run_num}.fasta"
        print(f"\nOutput filename: {final_output} (run {run_num})")

    # Print configuration
    print(f"\n{'='*70}")
    print(f"OPTIMIZED GENBANK FETCHER CONFIGURATION")
    print(f"{'='*70}")
    print(f"Input file: {args.input}")
    print(f"Output file: {final_output}")
    print(f"Total accessions: {len(accessions):,}")
    print(f"Batch size: {args.batch_size}")
    print(f"Base delay: {args.base_delay}s")
    print(f"Max delay: {args.max_delay}s")
    print(f"Max retries: {args.max_retries}")
    print(f"Email: {args.email}")
    if resume_mode:
        print(f"Resume file: {args.resume}")
    if args.compile_existing:
        print(f"Compile existing: Yes (from {args.output_dir or 'output/'})")
    print(f"{'='*70}")

    # Confirm with user for large datasets
    if len(accessions_to_fetch) > 10000 and not resume_mode and not args.compile_existing:
        print(f"\n⚠  Large dataset detected ({len(accessions_to_fetch):,} sequences)")
        print(f"Estimated time: ~{len(accessions_to_fetch) * 2 / 3600:.1f} hours (assuming 2 seconds per sequence)")
        response = input("Continue? (y/n): ")
        if response.lower() != 'y':
            print("Aborted by user")
            sys.exit(0)

    # Start retrieval
    start_time = time.time()
    
    # If compile_existing is used, we need to load existing sequences first
    if args.compile_existing and compiled_accessions:
        print(f"\n{'='*70}")
        print(f"LOADING EXISTING SEQUENCES")
        print(f"{'='*70}")
        existing_sequences = fetcher.load_existing_sequences(output_dir, compiled_accessions)
        print(f"Loaded {len(existing_sequences):,} existing sequences from FASTA files")
    
    success = fetcher.fetch_sequences(accessions_to_fetch, batch_size=args.batch_size, resume=resume_mode, 
                                     existing_fasta=args.resume, date_str=date_str)

    if success:
        # If compile_existing, combine existing + new sequences
        if args.compile_existing and compiled_accessions:
            fetcher.results = existing_sequences + fetcher.results
            print(f"\nCombined {len(existing_sequences):,} existing + {len(fetcher.results) - len(existing_sequences):,} new = {len(fetcher.results):,} total sequences")
        
        # Save results to final output file
        fetcher.save_results(final_output, resume=resume_mode)
        fetcher.display_sequence_stats()

        # Save failed accessions for retry
        if fetcher.failed_accessions:
            fetcher.save_failed_accessions()

        total_time = time.time() - start_time
        print(f"\nTotal execution time: {total_time:.1f} seconds ({total_time/60:.1f} minutes)")

        if fetcher.failed_accessions:
            print(f"\n⚠ Some sequences failed to retrieve ({len(set(fetcher.failed_accessions)):,})")
            print("  Run with --resume to retry failed sequences")
            sys.exit(1)
        else:
            print("\n✓ All sequences successfully retrieved!")
            sys.exit(0)
    else:
        print("\n✗ Failed to retrieve any sequences")
        sys.exit(1)

if __name__ == "__main__":
    import os  # Added for file size calculation
    main()
