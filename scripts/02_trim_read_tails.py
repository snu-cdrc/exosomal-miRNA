#!/usr/bin/env python3
"""Trim bases from the 3' end of FASTQ reads.

Removes the specified number of bases from both sequence and quality lines
in gzipped FASTQ files. Designed for removing NEXTflex random adapter
sequences from the 3' end of small RNA-seq reads.

Usage:
    python 02_trim_read_tails.py <input_dir> --trim-length 4
"""

import argparse
import gzip
import os
import sys


def trim_fastq(input_path, output_path, trim_length):
    """Trim `trim_length` bases from the 3' end of each read."""
    lines = []
    with gzip.open(input_path, "rb") as f:
        for i, line in enumerate(f):
            # Sequence (line 2) and quality (line 4) in each 4-line FASTQ record
            if i % 2 == 1:
                line = line[:-(trim_length + 1)] + b"\n"
            lines.append(line)

    with gzip.open(output_path, "wb") as f:
        for line in lines:
            f.write(line)

    print(f"  {os.path.basename(input_path)} -> {os.path.basename(output_path)}")


def main():
    parser = argparse.ArgumentParser(
        description="Trim bases from 3' end of FASTQ reads (NEXTflex random adapter removal)"
    )
    parser.add_argument("input_dir", help="Directory containing .fq.gz files")
    parser.add_argument(
        "--trim-length", type=int, required=True,
        help="Number of bases to trim from the 3' end"
    )
    args = parser.parse_args()

    if not os.path.isdir(args.input_dir):
        print(f"Error: {args.input_dir} is not a directory", file=sys.stderr)
        sys.exit(1)

    fq_files = sorted(f for f in os.listdir(args.input_dir) if f.endswith(".fq.gz"))
    if not fq_files:
        print(f"No .fq.gz files found in {args.input_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Trimming {args.trim_length}bp from 3' end of {len(fq_files)} files:")
    for fname in fq_files:
        input_path = os.path.join(args.input_dir, fname)
        output_path = os.path.join(args.input_dir, f"cut{args.trim_length}bp_{fname}")
        trim_fastq(input_path, output_path, args.trim_length)

    print("Done.")


if __name__ == "__main__":
    main()
