#!/usr/bin/env python3
"""Filter FASTA records for COMEBin and report simple record counts."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Iterator
from pathlib import Path
from typing import TextIO


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Write FASTA records longer than the COMEBin minimum length.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Input FASTA file.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Filtered output FASTA file.",
    )
    parser.add_argument(
        "--stats",
        type=Path,
        required=True,
        help="Output TSV with total, kept, and skipped counts.",
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=1000,
        help="Keep records strictly longer than this many bases.",
    )
    return parser.parse_args()


def iter_fasta(handle: TextIO) -> Iterator[tuple[str, str, str]]:
    """Yield FASTA records as header, identifier, and sequence text."""
    header: str | None = None
    seq_lines: list[str] = []

    for raw_line in handle:
        line = raw_line.rstrip("\n\r")
        if line.startswith(">"):
            if header is not None:
                seq = "".join(seq_lines)
                yield header, header[1:].split()[0], seq
            header = line
            seq_lines = []
            continue
        if header is None:
            raise ValueError("input does not look like FASTA: sequence before header")
        seq_lines.append(line.strip())

    if header is not None:
        seq = "".join(seq_lines)
        yield header, header[1:].split()[0], seq


def write_fasta_record(handle: TextIO, header: str, sequence: str) -> None:
    """Write one FASTA record with wrapped sequence lines."""
    handle.write(f"{header}\n")
    for start in range(0, len(sequence), 80):
        handle.write(f"{sequence[start:start + 80]}\n")


def filter_fasta(input_path: Path, output_path: Path, min_length: int) -> tuple[int, int]:
    """Filter FASTA records by length and return total and kept counts."""
    total = 0
    kept = 0
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with input_path.open("r", encoding="utf-8") as in_handle, output_path.open(
        "w", encoding="utf-8"
    ) as out_handle:
        for header, _record_id, sequence in iter_fasta(in_handle):
            total += 1
            if len(sequence) <= min_length:
                continue
            write_fasta_record(out_handle, header, sequence)
            kept += 1

    return total, kept


def write_stats(stats_path: Path, total: int, kept: int) -> None:
    """Write filtering counts as a small key-value TSV."""
    stats_path.parent.mkdir(parents=True, exist_ok=True)
    skipped = total - kept
    with stats_path.open("w", encoding="utf-8") as handle:
        handle.write(f"total\t{total}\n")
        handle.write(f"kept\t{kept}\n")
        handle.write(f"skipped\t{skipped}\n")


def main() -> int:
    """Run the COMEBin FASTA filter."""
    args = parse_args()
    if args.min_length < 0:
        print("--min-length must be non-negative", file=sys.stderr)
        return 1
    if not args.input.is_file():
        print(f"input FASTA not found: {args.input}", file=sys.stderr)
        return 1

    try:
        total, kept = filter_fasta(args.input, args.output, args.min_length)
        write_stats(args.stats, total, kept)
    except OSError as exc:
        print(f"FASTA filtering failed: {exc}", file=sys.stderr)
        return 1
    except ValueError as exc:
        print(f"FASTA filtering failed: {exc}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
