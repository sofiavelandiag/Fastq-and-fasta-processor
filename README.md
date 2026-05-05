# FASTX Preprocessor

A Python command-line tool that preprocesses DNA sequences from FASTA or FASTQ files. It supports three operations: reverse complement, hard trimming, and adaptor removal.

**Authors:** Sofía Velandia & Laura Stoma  
**Course:** MSc in Bioinformatics — Module I: Programming in Bioinformatics, October 2025

---

## Requirements

- Python 3.6 or higher
- No external libraries required (standard library only)

---

## Usage

```bash
python fastq_pp.py --input [FILE] --output [FILE] --operation [OPERATION]
```

Arguments can be provided in **any order**.

| Argument | Required | Description |
|---|---|---|
| `--input FILE` | Always | Input FASTA or FASTQ file |
| `--output FILE` | Always | Output file (must match input format) |
| `--operation` | Always | One of `rc`, `trim`, `adaptor-removal` |
| `--trim-left N` | Only for `trim` | Bases to remove from the left end |
| `--trim-right N` | Only for `trim` | Bases to remove from the right end |
| `--adaptor SEQ` | Only for `adaptor-removal` | Adaptor sequence to search for and remove |

---

## Supported Formats

The tool auto-detects the format from the file extension and the first non-empty line.

| Format | Valid extensions |
|---|---|
| FASTA | `.fasta`, `.fa`, `.FASTA`, `.FA` |
| FASTQ | `.fastq`, `.fq`, `.FASTQ`, `.FQ` |

> The input and output files must use the same format. Mixing `.fastq` input with `.fasta` output (or vice versa) will produce an error.

---

## Operations

### a) Reverse Complement (`--operation rc`)

Computes the reverse complement of every sequence. For FASTQ files, the quality string is reversed to match the new sequence orientation.

```bash
python fastq_pp.py --input test.fastq --output result.fastq --operation rc
```

**Example output:**
```
Detected input as FASTQ format.
Input and output are FASTQ files.
Reverse complement written to 'result.fastq'.

==================================================
Summary of Sequence Statistics — operation: rc
==================================================
Total reads processed: 100000
Total bases processed: 10000000
  A: 1500000 (15.00%)
  C: 4000000 (40.00%)
  G: 3500000 (35.00%)
  T: 1000000 (10.00%)
  N: 0 (0.00%)
```

---

### b) Trim (`--operation trim`)

Hard-trims a fixed number of bases from the left and/or right end of every sequence. At least one of `--trim-left` or `--trim-right` must be provided. The tool exits with an error if the combined trim length equals or exceeds the sequence length.

```bash
python fastq_pp.py --input test.fastq --output result.fastq --operation trim \
  --trim-left 20 --trim-right 10
```

**Example output:**
```
Detected input as FASTQ format.
Input and output are FASTQ files.
Trimmed sequences written to 'result.fastq'.

==================================================
Summary of Sequence Statistics — operation: trim
==================================================
Total reads processed: 100000
Total bases processed: 10000000
  A: 1500000 (15.00%)
  C: 4000000 (40.00%)
  G: 3500000 (35.00%)
  T: 1000000 (10.00%)
  N: 0 (0.00%)

Bases trimmed: 3000000
  A: 120000 (4.00%)
  C: 1500000 (50.00%)
  G: 1350000 (45.00%)
  T: 30000 (1.00%)
  N: 0 (0.00%)
```

---

### c) Adaptor Removal (`--operation adaptor-removal`)

Searches for the given adaptor sequence at the **start** of each read (case-insensitive). If found, it is trimmed off. Reads that do not start with the adaptor are kept unchanged.

```bash
python fastq_pp.py --input test.fastq --output result.fastq \
  --operation adaptor-removal \
  --adaptor TATAGA
```

**Example output:**
```
Detected input as FASTQ format.
Input and output are FASTQ files.
Sequences without adaptor written to 'result.fastq'.

==================================================
Summary of Sequence Statistics — operation: adaptor-removal
==================================================
Total reads processed: 100000
Total bases processed: 10000000
  A: 1500000 (15.00%)
  C: 4000000 (40.00%)
  G: 3500000 (35.00%)
  T: 1000000 (10.00%)
  N: 0 (0.00%)

Total adaptors found: 100
```

---

## Error Handling

| Error | Message |
|---|---|
| Missing required arguments | `Error: Missing required arguments.` |
| Unknown operation | `Error: Unknown operation. Use one of {rc, trim, adaptor-removal}` |
| `trim` with no trim values | `Error: For 'trim', provide --trim-left and/or --trim-right.` |
| `adaptor-removal` with no adaptor | `Error: For 'adaptor-removal', provide --adaptor SEQUENCE.` |
| Invalid input file extension | `Error: The input file '...' has a non-valid extension.` |
| Invalid output file extension | `Error: The output file '...' has a non-valid extension.` |
| Input/output format mismatch | `Error: The input and output files must have the same format.` |
| Trim removes entire sequence | `Error: Trim length exceeds or equals the read length.` |
| Adaptor longer than sequence | `Error: Adaptor length is equal to or longer than the sequence.` |

---

## Notes

- Multi-line FASTA sequences (sequence split across multiple lines) are fully supported.
- The format detection reads only the first non-empty line of the file, so it is fast regardless of file size.
- All operations stream through the file one record at a time, keeping memory usage low even for large files.
