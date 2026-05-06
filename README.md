# SeqRetrieve

A comprehensive toolkit for retrieving, sorting, aligning, filtering, and reporting GenBank sequences.

## Features

- **Sequence Retrieval**: Fetch sequences from NCBI GenBank with rate limiting and resume support
- **Metadata Retrieval**: Extract metadata (species, taxid, collection date, country, etc.)
- **Sorting**: Sort sequences and metadata by any field (taxid, organism, species, year, etc.)
- **Alignment**: Auto-align sequences grouped by taxonomic or temporal fields
- **Filtering**: Filter by sequence length, ambiguous base content (N count, IUPAC codes)
- **Reporting**: Generate HTML reports with classification statistics and filter results

## Installation

```bash
git clone https://github.com/Krysasp/SeqRetrieve.git
cd SeqRetrieve
pip install -r requirements.txt
pip install -e .  # optional: install as package
```

Dependencies:
- Python 3.7+
- Biopython
- numpy
- External alignment tools (optional): MUSCLE, Clustal Omega, or MAFFT

## Usage

The main CLI tool `SeqRetrieve` provides subcommands for each function:

### 1. Retrieve Sequences

```bash
./bin/SeqRetrieve retrieve -i accessions.csv -o sequences.fasta -e your@email.com
```

Options:
- `-i, --input`: Input CSV file with accession numbers
- `-o, --output`: Output FASTA file
- `-e, --email`: Your email (required by NCBI)
- `-r, --resume`: Resume from existing FASTA file
- `--batch-size`: Batch size (default 50)
- `--max-retries`: Maximum retries per batch (default 5)

### 2. Retrieve Metadata

```bash
./bin/SeqRetrieve metadata -i accessions.csv -o metadata.tsv -e your@email.com
```

Or extract accessions from existing FASTA:
```bash
./bin/SeqRetrieve metadata -f sequences.fasta -o metadata.tsv -e your@email.com
```

### 3. Sort Sequences and Metadata

```bash
./bin/SeqRetrieve sort -f sequences.fasta -m metadata.tsv --field species \
  -of sorted.fasta -om sorted.tsv
```

Supported fields: `taxid`, `organism`, `species`, `genus`, `collection_date`, `country`, `genotype`, `host`, `length`

### 4. Align Sequences

```bash
./bin/SeqRetrieve align -f sequences.fasta -m metadata.tsv --field species -o aligned/
```

Requires one of: `muscle`, `clustalo`, or `mafft` in PATH.

### 5. Filter Sequences

```bash
./bin/SeqRetrieve filter -f sequences.fasta -o output/ \
  --min-length-pct 80 --max-ambig-pct 5
```

Options:
- `--min-length-pct`: Minimum length as % of longest sequence
- `--max-ambig-count`: Maximum count of ambiguous bases (N + IUPAC)
- `--max-ambig-pct`: Maximum percentage of ambiguous bases
- `--filtered-subdir`: Subdirectory for filtered sequences (default: `filtered`)

### 6. Generate HTML Report

```bash
./bin/SeqRetrieve report -o report.html -m metadata.tsv \
  --filter-report output/filter_report.tsv
```

## Example Workflow

```bash
# 1. Retrieve sequences
./bin/SeqRetrieve retrieve -i accessions.csv -o seqs.fasta -e youremail@example.com

# 2. Retrieve metadata
./bin/SeqRetrieve metadata -f seqs.fasta -o meta.tsv -e youremail@example.com

# 3. Filter sequences (remove short/ambiguous)
./bin/SeqRetrieve filter -f seqs.fasta -o filtered/ --min-length-pct 80 --max-ambig-pct 5

# 4. Sort by species
./bin/SeqRetrieve sort -f filtered/seqs.fasta -m meta.tsv --field species \
  -of sorted/seqs.fasta -om sorted/meta.tsv

# 5. Align by species
./bin/SeqRetrieve align -f sorted/seqs.fasta -m sorted/meta.tsv \
  --field species -o aligned/

# 6. Generate report
./bin/SeqRetrieve report -o report.html -m sorted/meta.tsv \
  --filter-report filtered/filter_report.tsv
```

## Resume Support

Both `retrieve` and `metadata` commands support resume mode to continue interrupted runs:

```bash
./bin/SeqRetrieve retrieve -i accessions.csv -o seqs.fasta -e youremail@example.com \
  -r seqs.fasta
```

## License

MIT License
