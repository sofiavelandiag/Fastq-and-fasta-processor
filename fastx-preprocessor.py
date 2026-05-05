#!/usr/bin/env python3

# -----------------------------------------------------------------------------
# DNA Sequence Processor (FASTA/FASTQ)
# -----------------------------------------------------------------------------
# This program performs 3 main operations on DNA sequences:
#   1) Reverse complement                        (--operation rc)
#   2) Trimming of bases from one or both ends   (--operation trim)
#   3) Adaptor sequence removal (left end only)  (--operation adaptor-removal)
#
# It reads sequences from FASTA or FASTQ files (automatically detected)
# and writes summary statistics (counts and base percentages).
# -----------------------------------------------------------------------------

import sys

# --------------------------------------------------------------
# ARGUMENT PARSING
# --------------------------------------------------------------
# Default values for all expected command-line arguments.
# trim_left and trim_right hold numerical values; others are strings.
input_file  = ""
output_file = ""
operation   = ""
trim_left   = 0
trim_right  = 0
adaptor_seq = ""

# Reads command-line arguments provided by the user.
# Assigns values to the variables that follow each flag.
for i in range(1, len(sys.argv)):
    arg = sys.argv[i]
    if arg == "--input" and i + 1 < len(sys.argv):
        input_file  = sys.argv[i + 1]
    elif arg == "--output" and i + 1 < len(sys.argv):
        output_file = sys.argv[i + 1]
    elif arg == "--operation" and i + 1 < len(sys.argv):
        operation   = sys.argv[i + 1]
    elif arg == "--trim-left" and i + 1 < len(sys.argv):
        trim_left   = int(sys.argv[i + 1])
    elif arg == "--trim-right" and i + 1 < len(sys.argv):
        trim_right  = int(sys.argv[i + 1])
    elif arg == "--adaptor" and i + 1 < len(sys.argv):
        adaptor_seq = sys.argv[i + 1]


# --------------------------------------------------------------
# ARGUMENT VALIDATION
# --------------------------------------------------------------
# Ensures all required options are provided correctly before running.

# Check that --input, --output and --operation were all supplied
if not input_file or not output_file or not operation:
    print("Error: Missing required arguments.")
    print("Usage: --input FILE --output FILE --operation {rc,trim,adaptor-removal}")
    sys.exit(1)

# Ensure operation type is valid
if operation not in {"rc", "trim", "adaptor-removal"}:
    print("Error: Unknown operation. Use one of {rc, trim, adaptor-removal}")
    sys.exit(1)

# Ensure trim operation has at least one side specified
if operation == "trim" and (trim_left == 0 and trim_right == 0):
    print("Error: For 'trim', provide --trim-left and/or --trim-right.")
    sys.exit(1)

# Ensure adaptor removal has a sequence provided
if operation == "adaptor-removal" and adaptor_seq == "":
    print("Error: For 'adaptor-removal', provide --adaptor SEQUENCE.")
    sys.exit(1)


# --------------------------------------------------------------
# FORMAT DETECTION & VALIDATION
# --------------------------------------------------------------

def detect_format(path: str) -> str:
    """
    Detect whether a file is in FASTA or FASTQ format.

    Checks the file extension first, then inspects the first
    non-empty line: '>' indicates FASTA and '@' indicates FASTQ.

    Parameters:
        path (str): Path to the input file.

    Returns:
        str: 'fasta' or 'fastq'.

    Exits:
        sys.exit(1) if the extension is invalid or the format
        cannot be determined from the file content.
    """
    valid_extensions = (".fasta", ".fastq", ".fa", ".fq",
                        ".FASTA", ".FASTQ", ".FA", ".FQ")

    if not path.endswith(valid_extensions):
        print(f"Error: The input file '{path}' has a non-valid extension.")
        print(f"Allowed extensions are: {', '.join(valid_extensions)}")
        sys.exit(1)

    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("@"):
                print("Detected input as FASTQ format.")
                return "fastq"
            if s.startswith(">"):
                print("Detected input as FASTA format.")
                return "fasta"
            break  # First non-empty line was neither '@' nor '>'

    print("Error: Could not detect FASTA/FASTQ "
          "(first non-empty line does not start with '@' or '>').")
    sys.exit(1)


def validate_output_format(path: str):
    """
    Validate that the output file path has a recognised sequence file extension.

    Parameters:
        path (str): Path to the output file.

    Exits:
        sys.exit(1) if the extension is not recognised.
    """
    valid_extensions = (".fasta", ".fastq", ".fa", ".fq",
                        ".FASTA", ".FASTQ", ".FA", ".FQ")
    if not path.endswith(valid_extensions):
        print(f"Error: The output file '{path}' has a non-valid extension.")
        print(f"Allowed extensions are: {', '.join(valid_extensions)}")
        sys.exit(1)


def validate_same_format(input_path: str, output_path: str):
    """
    Ensure that the input and output files share the same sequence format
    (both FASTA or both FASTQ).

    Parameters:
        input_path  (str): Path to the input file.
        output_path (str): Path to the output file.

    Exits:
        sys.exit(1) if the two files belong to different formats.
    """
    fastq_ext = (".fastq", ".FASTQ", ".fq", ".FQ")
    fasta_ext = (".fasta", ".FASTA", ".fa", ".FA")

    if input_path.endswith(fastq_ext) and output_path.endswith(fastq_ext):
        print("Input and output are FASTQ files.")
    elif input_path.endswith(fasta_ext) and output_path.endswith(fasta_ext):
        print("Input and output are FASTA files.")
    else:
        print("Error: The input and output files must have the same format.")
        sys.exit(1)


# Run format checks
fmt = detect_format(input_file)
validate_output_format(output_file)
validate_same_format(input_file, output_file)


# --------------------------------------------------------------
# SEQUENCE STATISTICS
# --------------------------------------------------------------

class SequenceStats:
    """
    Accumulate and report nucleotide statistics across processed sequences.

    Attributes:
        reads (int): Number of sequences processed.
        bases (int): Total number of bases counted.
        A (int): Count of adenine bases.
        C (int): Count of cytosine bases.
        G (int): Count of guanine bases.
        T (int): Count of thymine bases.
        N (int): Count of ambiguous bases.
    """

    def __init__(self):
        """Initialise all counters to zero."""
        self.reads = 0
        self.bases = 0
        self.A = 0
        self.C = 0
        self.G = 0
        self.T = 0
        self.N = 0

    def count_bases(self, seq: str):
        """
        Count each nucleotide in a sequence and update cumulative totals.

        Parameters:
            seq (str): A DNA sequence string (upper or lowercase).
        """
        self.A += seq.count("A") + seq.count("a")
        self.C += seq.count("C") + seq.count("c")
        self.G += seq.count("G") + seq.count("g")
        self.T += seq.count("T") + seq.count("t")
        self.N += seq.count("N") + seq.count("n")
        self.bases += len(seq)
        self.reads += 1

    def pct(self, x: int) -> str:
        """
        Calculate the percentage of a base count relative to total bases.

        Parameters:
            x (int): Count of a specific base.

        Returns:
            str: Percentage formatted to two decimal places, or '0.00'
                 if no bases have been counted yet.
        """
        return f"{(100.0 * x / self.bases):.2f}" if self.bases else "0.00"

    def print_stats(self):
        """Print a formatted summary of accumulated sequence statistics."""
        print("\n" + "=" * 50)
        print(f"Summary of Sequence Statistics — operation: {operation}")
        print("=" * 50)
        print(f"Total reads processed: {self.reads}")
        print(f"Total bases processed: {self.bases}")
        print(f"  A: {self.A} ({self.pct(self.A)}%)")
        print(f"  C: {self.C} ({self.pct(self.C)}%)")
        print(f"  G: {self.G} ({self.pct(self.G)}%)")
        print(f"  T: {self.T} ({self.pct(self.T)}%)")
        print(f"  N: {self.N} ({self.pct(self.N)}%)")


# --------------------------------------------------------------
# FILE ITERATOR
# --------------------------------------------------------------

def iter_file(path: str, format_type: str):
    """
    Generator that yields records from a FASTA or FASTQ file.

    For FASTQ files, each call yields a 4-tuple:
        (header, sequence, plus_line, quality)
    For FASTA files, each call yields a 2-tuple:
        (header, sequence)

    All strings are returned with trailing newlines preserved for
    FASTQ (so the caller can strip selectively) and stripped for FASTA.

    Parameters:
        path        (str): Path to the input file.
        format_type (str): 'fastq' or 'fasta'.
    """
    with open(path, "r") as f:

        # ── FASTQ: strict 4-line records ──────────────────────────
        if format_type == "fastq":
            while True:
                h = f.readline()
                if not h:       # End of file
                    return
                s = f.readline()
                p = f.readline()
                q = f.readline()
                yield h, s, p, q

        # ── FASTA: header + one-or-more sequence lines ────────────
        elif format_type == "fasta":
            header    = None
            seq_parts = []

            for line in f:
                if line.startswith(">"):
                    # Yield the previous record before starting a new one
                    if header is not None:
                        yield header, "".join(seq_parts)
                    # Start a new record
                    header    = line.rstrip("\n")
                    seq_parts = []
                else:
                    # Accumulate sequence lines (handles multi-line FASTA)
                    seq_parts.append(line.strip())

            # Yield the final record (would otherwise be silently dropped)
            if header is not None:
                yield header, "".join(seq_parts)


# --------------------------------------------------------------
# OPERATION 1: REVERSE COMPLEMENT
# --------------------------------------------------------------

# Translation table: each base maps to its Watson-Crick complement
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s: str) -> str:
    """
    Return the reverse complement of a DNA sequence.

    Parameters:
        s (str): A DNA sequence string.

    Returns:
        str: The reverse complement (bases complemented, order reversed).
    """
    return s.translate(COMP)[::-1]


def run_rc(input_file: str, output_file: str, fmt: str):
    """
    Compute the reverse complement of every sequence in a FASTA/FASTQ file.

    For FASTQ files the quality string is reversed to match the new sequence
    orientation (quality scores cannot be complemented).

    Parameters:
        input_file  (str): Path to the input FASTA or FASTQ file.
        output_file (str): Path to the output file.
        fmt         (str): Detected format — 'fasta' or 'fastq'.
    """
    stats = SequenceStats()

    with open(output_file, "w") as fo:

        if fmt == "fastq":
            for h, s, p, q in iter_file(input_file, fmt):
                s = s.strip("\n")
                q = q.strip("\n")

                rc_seq  = revcomp(s)
                rc_qual = q[::-1]   # Reverse quality to match flipped sequence

                fo.write(h)
                fo.write(rc_seq  + "\n")
                fo.write(p)
                fo.write(rc_qual + "\n")

                stats.count_bases(rc_seq)

        elif fmt == "fasta":
            for header, seq in iter_file(input_file, fmt):
                seq    = seq.strip("\n")
                rc_seq = revcomp(seq)

                fo.write(header + "\n")
                fo.write(rc_seq + "\n")

                stats.count_bases(rc_seq)

    print(f"Reverse complement written to '{output_file}'.")
    stats.print_stats()


# --------------------------------------------------------------
# OPERATION 2: TRIM BASES
# --------------------------------------------------------------

def run_trim(input_file: str, output_file: str, fmt: str,
             trim_left: int, trim_right: int):
    """
    Hard-trim a fixed number of bases from the left and/or right of every
    sequence in a FASTA or FASTQ file.

    The function exits with an error if the combined trim length is greater
    than or equal to the sequence length, which would produce an empty read.

    Parameters:
        input_file  (str): Path to the input FASTA or FASTQ file.
        output_file (str): Path to the output file.
        fmt         (str): Detected format — 'fasta' or 'fastq'.
        trim_left   (int): Number of bases to remove from the left end.
        trim_right  (int): Number of bases to remove from the right end.
    """
    stats        = SequenceStats()
    trim_stats   = SequenceStats()   # Tracks the bases that were removed

    with open(output_file, "w") as fo:

        if fmt == "fastq":
            for h, s, p, q in iter_file(input_file, fmt):
                s = s.strip("\n")
                q = q.strip("\n")
                L = len(s)

                if (trim_left + trim_right) >= L:
                    print("Error: Trim length exceeds or equals the read length.")
                    sys.exit(1)

                end    = L - trim_right if trim_right > 0 else None
                kept   = s[trim_left:end]
                kept_q = q[trim_left:end]

                # Track trimmed bases for the summary
                left_trim  = s[:trim_left]
                right_trim = s[L - trim_right:] if trim_right > 0 else ""
                if left_trim:
                    trim_stats.count_bases(left_trim)
                if right_trim:
                    trim_stats.count_bases(right_trim)

                fo.write(h)
                fo.write(kept   + "\n")
                fo.write(p)
                fo.write(kept_q + "\n")

                stats.count_bases(kept)

        elif fmt == "fasta":
            for header, s in iter_file(input_file, fmt):
                s = s.strip("\n")
                L = len(s)

                if (trim_left + trim_right) >= L:
                    print("Error: Trim length exceeds or equals the read length.")
                    sys.exit(1)

                end  = L - trim_right if trim_right > 0 else None
                kept = s[trim_left:end]

                left_trim  = s[:trim_left]
                right_trim = s[L - trim_right:] if trim_right > 0 else ""
                if left_trim:
                    trim_stats.count_bases(left_trim)
                if right_trim:
                    trim_stats.count_bases(right_trim)

                fo.write(header + "\n")
                fo.write(kept   + "\n")

                stats.count_bases(kept)

    print(f"Trimmed sequences written to '{output_file}'.")
    stats.print_stats()

    # Print trimmed-bases summary (mirrors the assignment's expected output)
    print(f"\nBases trimmed: {trim_stats.bases}")
    if trim_stats.bases:
        print(f"  A: {trim_stats.A} ({trim_stats.pct(trim_stats.A)}%)")
        print(f"  C: {trim_stats.C} ({trim_stats.pct(trim_stats.C)}%)")
        print(f"  G: {trim_stats.G} ({trim_stats.pct(trim_stats.G)}%)")
        print(f"  T: {trim_stats.T} ({trim_stats.pct(trim_stats.T)}%)")
        print(f"  N: {trim_stats.N} ({trim_stats.pct(trim_stats.N)}%)")


# --------------------------------------------------------------
# OPERATION 3: ADAPTOR REMOVAL
# --------------------------------------------------------------

def run_adaptor(input_file: str, output_file: str, fmt: str, adaptor_seq: str):
    """
    Search for an adaptor sequence at the start of each read and remove it.

    The match is case-insensitive. Reads that do not begin with the adaptor
    are written to the output unchanged. The function exits with an error if
    the adaptor is longer than or equal to a read.

    Parameters:
        input_file  (str): Path to the input FASTA or FASTQ file.
        output_file (str): Path to the output file.
        fmt         (str): Detected format — 'fasta' or 'fastq'.
        adaptor_seq (str): The adaptor sequence to search for and remove.
    """
    stats          = SequenceStats()
    ADAP           = adaptor_seq.upper()    # Normalise adaptor to uppercase
    adaptors_found = 0

    with open(output_file, "w") as fo:

        if fmt == "fastq":
            for h, s, p, q in iter_file(input_file, fmt):
                s = s.strip("\n")
                q = q.strip("\n")

                if len(ADAP) >= len(s):
                    print("Error: Adaptor length is equal to or longer than the sequence.")
                    sys.exit(1)

                if s.upper().startswith(ADAP):
                    adaptors_found += 1
                    trimmed   = s[len(ADAP):]
                    trimmed_q = q[len(ADAP):]
                else:
                    trimmed   = s
                    trimmed_q = q

                fo.write(h)
                fo.write(trimmed   + "\n")
                fo.write(p)
                fo.write(trimmed_q + "\n")

                stats.count_bases(trimmed)

        elif fmt == "fasta":
            for header, s in iter_file(input_file, fmt):
                s = s.strip("\n")

                if len(ADAP) >= len(s):
                    print("Error: Adaptor length is equal to or longer than the sequence.")
                    sys.exit(1)

                if s.upper().startswith(ADAP):
                    adaptors_found += 1
                    trimmed = s[len(ADAP):]
                else:
                    trimmed = s

                fo.write(header  + "\n")
                fo.write(trimmed + "\n")

                stats.count_bases(trimmed)

    print(f"Sequences without adaptor written to '{output_file}'.")
    stats.print_stats()
    print(f"\nTotal adaptors found: {adaptors_found}")


# --------------------------------------------------------------
# MAIN OPERATION DISPATCH
# --------------------------------------------------------------

# Direct the script to the appropriate function based on the user's request.
if operation == "rc":
    run_rc(input_file, output_file, fmt)

elif operation == "trim":
    run_trim(input_file, output_file, fmt, trim_left, trim_right)

elif operation == "adaptor-removal":
    run_adaptor(input_file, output_file, fmt, adaptor_seq)