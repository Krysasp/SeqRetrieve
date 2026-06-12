# SeqRetrieve

A comprehensive toolkit for retrieving, sorting, aligning, filtering, and reporting GenBank sequences.

## Features

- **Sequence Retrieval**: Fetch sequences from NCBI GenBank with batching, retry logic, and resume support
- **Metadata Retrieval**: Extract metadata (species, taxid, collection date, country, host, genotype, etc.)
- **Sorting**: Sort sequences and metadata by any field (taxid, organism, species, year, country, host, etc.)
- **Alignment**: Auto-align sequences grouped by taxonomic or temporal fields using MUSCLE, Clustal Omega, or MAFFT
- **Filtering**: Filter by sequence length, ambiguous base content (N count, IUPAC codes)
- **Reporting**: Generate HTML reports with classification statistics and filter results
- **Compile Existing FASTA**: Combine accession IDs from multiple existing FASTA files to avoid redundant downloads

## Installation

```bash
git clone https://github.com/Krysasp/SeqRetrieve.git
cd SeqRetrieve
pip install -r requirements.txt
pip install -e .  # optional: install as package
```

### Dependencies

- Python 3.7+
- Biopython >= 1.80
- numpy >= 1.20
- External alignment tools (optional): MUSCLE, Clustal Omega, or MAFFT

### Project Structure

```
SeqRetrieve/
├── bin/
│   └── retrieveseq          # Main CLI executable
├── src/
│   ├── retrieve_seq.py      # Sequence retrieval module
│   ├── retrieve_metadata.py # Metadata extraction module
│   ├── sort_seqs.py         # Sorting module
│   ├── align_seqs.py        # Alignment module
│   ├── filter_seqs.py       # Filtering module
│   └── generate_report.py   # HTML report generation module
├── tests/                   # Test files
├── docs/                    # Documentation
├── requirements.txt         # Python dependencies
├── setup.py                 # Package installation
├── .gitignore              # Git ignore patterns
└── LICENSE                 # MIT License
```

## Usage

The main CLI tool `SeqRetrieve` provides subcommands for each function:

```bash
./bin/retrieveseq COMMAND [OPTIONS]
```

### Available Commands

| Command | Description |
|---------|-------------|
| `retrieve` | Download sequences from GenBank by accession numbers |
| `metadata` | Extract metadata from GenBank records |
| `sort` | Sort FASTA sequences and metadata TSV by any field |
| `align` | Align sequences grouped by field using MUSCLE/Clustal Omega/MAFFT |
| `filter` | Filter sequences by length and ambiguous base content |
| `report` | Generate HTML report with statistics and classification summaries |

### 1. Retrieve Sequences

```bash
./bin/retrieveseq retrieve -i accessions.csv -o sequences.fasta -e your@email.com
```

**Options:**
- `-i, --input`: Input CSV file with accession numbers (required)
- `-o, --output`: Output FASTA file (required, uses date-based naming: YYYY-MM-DD_name_runN.fasta)
- `-e, --email`: Email address for NCBI Entrez API (required)
- `-c, --column`: Column name containing accession IDs (default: first column)
- `-r, --resume`: Resume from existing FASTA file (skips already downloaded accessions)
- `--batch-size`: Batch size (default: 50, max recommended: 50)
- `--max-retries`: Maximum retry attempts per batch (default: 5)
- `--compile-existing`: Compile accession IDs from existing FASTA files
- `--output-dir`: Directory containing existing FASTA files (default: output directory)

### 2. Retrieve Metadata

```bash
./bin/retrieveseq metadata -i accessions.csv -o metadata.tsv -e your@email.com
```

Or extract accessions from existing FASTA:
```bash
./bin/retrieveseq metadata -f sequences.fasta -o metadata.tsv -e your@email.com
```

**Options:**
- `-i, --input`: Input CSV file with accession numbers
- `-f, --fasta`: Input FASTA file to extract accessions from
- `-o, --output`: Output TSV file for metadata (required)
- `-e, --email`: Email address for NCBI Entrez API (required)
- `-c, --column`: Column name containing accession IDs (default: first column)
- `-r, --resume`: Resume from existing metadata TSV
- `--batch-size`: Batch size (default: 50)
- `--max-retries`: Maximum retry attempts per batch (default: 5)

**Extracted metadata fields:** accession, version, length, taxid, organism, species, genus, family, collection_date, country, isolate, strain, host, genotype, submitter, seq_tech

### 3. Sort Sequences and Metadata

```bash
./bin/retrieveseq sort -f sequences.fasta -m metadata.tsv --field species \
  -of sorted.fasta -om sorted.tsv
```

**Options:**
- `-f, --fasta`: Input FASTA file (required)
- `-m, --metadata`: Input metadata TSV file (required)
- `--field`: Field to sort by (required)
- `-of, --output-fasta`: Output sorted FASTA file (required)
- `-om, --output-metadata`: Output sorted metadata TSV file (required)

**Supported sort fields:** `taxid`, `organism`, `species`, `genus`, `family`, `country`, `host`, `genotype`, `collection_date`, `length`

### 4. Align Sequences

```bash
./bin/retrieveseq align -f sequences.fasta -m metadata.tsv --field species -o aligned/
```

**Options:**
- `-f, --fasta`: Input FASTA file (required)
- `-m, --metadata`: Input metadata TSV file (required)
- `--field`: Field to group by (required)
- `-o, --output-dir`: Output directory for aligned FASTA files (required)
- `--tool`: Alignment tool: `muscle`, `clustalo`, or `mafft` (auto-detected if not specified)
- `--combined-output`: Optional combined aligned FASTA file

### 5. Filter Sequences

```bash
./bin/retrieveseq filter -f sequences.fasta -o output/ \
  --min-length-pct 80 --max-ambig-pct 5
```

**Options:**
- `-f, --fasta`: Input FASTA file (required)
- `-o, --output-dir`: Output directory for passed sequences (required)
- `--filtered-subdir`: Subdirectory for filtered sequences (default: `filtered`)
- `--min-length-pct`: Minimum length as % of longest sequence
- `--max-ambig-count`: Maximum count of ambiguous bases (N + IUPAC)
- `--max-ambig-pct`: Maximum percentage of ambiguous bases

### 6. Generate HTML Report

```bash
./bin/retrieveseq report -o report.html -m metadata.tsv \
  --filter-report output/filter_report.tsv
```

**Options:**
- `-o, --output`: Output HTML report file (required)
- `-f, --fasta`: Input FASTA file (for sequence statistics)
- `-m, --metadata`: Input metadata TSV file (for classification summaries)
- `--filter-report`: Filter report TSV file (for filtering details)
- `--filter-stats`: JSON file with filter statistics (optional)

## Example Workflow

Complete workflow for retrieving, filtering, sorting, aligning, and reporting on GenBank sequences:

```bash
# 1. Retrieve sequences from GenBank
./bin/retrieveseq retrieve -i accessions.csv -o seqs.fasta -e youremail@example.com

# 2. Retrieve metadata
./bin/retrieveseq metadata -f seqs.fasta -o meta.tsv -e youremail@example.com

# 3. Filter sequences (remove short/ambiguous)
./bin/retrieveseq filter -f seqs.fasta -o filtered/ \
  --min-length-pct 80 --max-ambig-pct 5

# 4. Sort by species
./bin/retrieveseq sort -f filtered/seqs.fasta -m meta.tsv --field species \
  -of sorted/seqs.fasta -om sorted/meta.tsv

# 5. Align by species
./bin/retrieveseq align -f sorted/seqs.fasta -m sorted/meta.tsv \
  --field species -o aligned/

# 6. Generate report
./bin/retrieveseq report -o report.html -m sorted/meta.tsv \
  --filter-report filtered/filter_report.tsv
```

## Advanced Features

### Resume Support

Both `retrieve` and `metadata` commands support resume mode to continue interrupted runs:

```bash
./bin/retrieveseq retrieve -i accessions.csv -o seqs.fasta -e youremail@example.com \
  -r seqs.fasta
```

This skips accession IDs that are already present in the output file.

### Compile Existing FASTA Files

The `--compile-existing` option allows you to combine accession IDs from multiple existing FASTA files with new accessions from a CSV file. This prevents redundant downloads when working with partial retrievals.

```bash
./bin/retrieveseq retrieve -i accessions.csv -o final_sequences.fasta \
  -e youremail@example.com --compile-existing --output-dir output/
```

**Workflow:**
1. Scans the specified `--output-dir` for all FASTA files
2. Extracts all accession IDs from existing files
3. Compares with CSV input to identify missing accessions
4. Downloads only the missing sequences
5. Creates consolidated output with date-based naming: `YYYY-MM-DD_<original_filename>.fasta`

**Example:**
```bash
# Multiple partial retrievals
./bin/retrieveseq retrieve -i batch1.csv -o output/batch1.fasta -e youremail@example.com
./bin/retrieveseq retrieve -i batch2.csv -o output/batch2.fasta -e youremail@example.com
./bin/retrieveseq retrieve -i batch3.csv -o output/batch3.fasta -e youremail@example.com

# Compile all and download only missing sequences
./bin/retrieveseq retrieve -i accessions.csv -o final_sequences.fasta \
  -e youremail@example.com --compile-existing --output-dir output/
```

## Quick Start

### Common Use Cases

**Retrieve a batch of sequences:**
```bash
./bin/retrieveseq retrieve -i my_accessions.csv -o my_sequences.fasta -e myemail@example.com
```

**Extract metadata from existing FASTA:**
```bash
./bin/retrieveseq metadata -f sequences.fasta -o metadata.tsv -e myemail@example.com
```

**Filter low-quality sequences:**
```bash
./bin/retrieveseq filter -f sequences.fasta -o clean/ \
  --min-length-pct 90 --max-ambig-pct 2
```

**Align sequences by taxonomy:**
```bash
./bin/retrieveseq align -f sequences.fasta -m metadata.tsv \
  --field taxid -o alignments/
```

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Run tests
5. Submit a pull request

## License

MIT License - See [LICENSE](LICENSE) file for details

## Author

SeqRetrieve toolkit developed by RetrieveSeq Contributors

## Acknowledgments

- NCBI for providing the Entrez API
- Biopython developers for sequence parsing utilities
- Alignment tool developers: MUSCLE, Clustal Omega, MAFFT
