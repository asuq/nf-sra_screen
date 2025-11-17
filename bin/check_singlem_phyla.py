#!/usr/bin/env python3
"""
Check presence of target phyla (GTDB style) in a single-sample taxonomy TSV.

Exit status:
0 : at least one target phylum is present (i.e., some row has present == 'yes')
1 : internal error (I/O problems, malformed input, unexpected exceptions)
2 : no target phylum is present (i.e., all rows have present == 'no' or no rows)
Usage:
    python check_phyla.py \
        --input sample_taxonomy.tsv \
        --phyla-list phyla_to_check.txt \
        --output phyla_presence.tsv \
        --verbose

CLI:
    -i, --input : TSV output from singlem pipe (header: 'sample\tcoverage\ttaxonomy')
        Example taxonomy cell: 'Root; d__Bacteria; p__Bacillota; c__Bacilli; ...'

    -p, --phyla-list: Phyla list text file (one per line) with a header line 'phyla'.
        Entries are GTDB-style (e.g. 'p__Bacillota'). Matching is case-insensitive
        with normalisation that strips the 'p__' prefix on both sides.

    -o, --output: TSV with columns (phyla, present)
        where 'present' is 'yes' or 'no' for each requested phylum

Required Packages: pandas

Author: Akito Shima (asuq)
Email: asuq.4096@gmail.com
"""

import argparse
import logging
import sys
from collections.abc import Iterable, Sequence
from pathlib import Path

import pandas as pd


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse CLI arguments."""
    p = argparse.ArgumentParser(
        description=(
            "Check whether a single-sample taxonomy TSV contains target phyla "
            "(GTDB style). Matching is case-insensitive with 'p__' normalisation."
        )
    )
    p.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Path to the input sample TSV with a 'taxonomy' column.",
    )
    p.add_argument(
        "-p",
        "--phyla-list",
        type=Path,
        required=True,
        help="Path to phyla list text file (header 'phyla', one per line).",
    )
    p.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Path to output TSV (columns: phyla, present).",
    )
    p.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose (DEBUG) logging.",
    )
    return p.parse_args(argv)


def setup_logging(verbose: bool) -> None:
    """Initialise logging with INFO or DEBUG level."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(levelname)s: %(message)s")


def normalise_phylum_label(label: str) -> str:
    """
    Normalise a GTDB-style phylum label for matching.

    Behaviour:
    - Strips surrounding whitespace.
    - Removes a leading 'p__' (case-insensitive) if present.
    - Lower-cases the remainder.

    Examples
    --------
    'p__Bacillota' -> 'bacillota'
    'P__PROTEOBACTERIA' -> 'proteobacteria'
    'bacillota' -> 'bacillota'
    """
    s = label.strip()
    if s.lower().startswith("p__"):
        s = s[3:]
    return s.strip().lower()


def extract_phylum_from_taxonomy(taxonomy: str) -> str | None:
    """
    Extract the phylum (normalised; no 'p__', lower-case) from a taxonomy path.

    The taxonomy string is expected to be semicolon-separated. We scan tokens
    to find the first one that starts with 'p__' (case-insensitive). If found,
    return the normalised phylum; otherwise return None.
    """
    if not isinstance(taxonomy, str):
        return None

    # Tokenise on ';' and trim whitespace; robust to extra spaces.
    for token in (t.strip() for t in taxonomy.split(";")):
        if token.lower().startswith("p__"):
            return normalise_phylum_label(token)
    return None


def read_observed_phyla_from_input(input_path: Path) -> set[str]:
    """
    Read the input TSV and return the set of observed phyla (normalised).

    Only the 'taxonomy' column is used. Rows without a phylum are ignored.
    """
    logging.info("Reading taxonomy TSV: %s", input_path)
    df = pd.read_csv(input_path, sep="\t", dtype=str, usecols=["taxonomy"])
    observed: set[str] = set()
    for tax in df["taxonomy"].dropna():
        ph = extract_phylum_from_taxonomy(tax)
        if ph:
            observed.add(ph)
    logging.info("Observed %d unique phyla in file.", len(observed))
    logging.debug("Observed phyla (normalised): %s", sorted(observed))
    return observed


def read_target_phyla(list_path: Path) -> list[tuple[str, str]]:
    """
    Read the phyla list file and return pairs of (original_label, normalised_key).

    The file contains a single column with header 'phyla'. Lines may include a
    'p__' prefix (GTDB style). Empty lines are ignored. Order is preserved.
    """
    logging.info("Reading phyla list: %s", list_path)
    raw_text = list_path.read_text(encoding="utf-8")

    targets: list[tuple[str, str]] = []
    for i, raw in enumerate(raw_text.splitlines()):
        line = raw.strip().lstrip("\ufeff")  # handle any BOM on first line
        if not line:
            continue
        if i == 0 and line.lower() == "phyla":
            continue
        normalised = normalise_phylum_label(line)
        targets.append((line, normalised))

    if not targets:
        logging.warning("No phyla found in list file after processing header.")
    else:
        logging.info("Loaded %d target phyla.", len(targets))
        logging.debug("Targets (original -> normalised): %s", targets)
    return targets


def compute_presence_records(
    targets: Sequence[tuple[str, str]], observed: set[str]
) -> list[tuple[str, str]]:
    """
    For each target phylum, compute presence ('yes'/'no') against the observed set.

    Returns a list of (original_label, present_flag).
    """
    records: list[tuple[str, str]] = []
    for original, key in targets:
        present = "yes" if key in observed else "no"
        records.append((original, present))
    return records


def write_output(records: Iterable[tuple[str, str]], output_path: Path) -> None:
    """Write the results to a TSV with columns 'phyla' and 'present'."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df_out = pd.DataFrame(records, columns=["phyla", "present"])
    df_out.to_csv(output_path, sep="\t", index=False)
    logging.info("Wrote results: %s", output_path)


def main(argv: Sequence[str] | None = None) -> None:
    """CLI entry point."""
    try:
        args = parse_args(argv)
    except SystemExit as se:
        if se.code != 0:
            sys.exit(1)
        raise

    setup_logging(args.verbose)

    try:
        if not args.input.exists():
            raise FileNotFoundError(f"Input TSV not found: {args.input}")
        if not args.phyla_list.exists():
            raise FileNotFoundError(f"Phyla list not found: {args.phyla_list}")

        observed = read_observed_phyla_from_input(args.input)
        targets = read_target_phyla(args.phyla_list)
        records = compute_presence_records(targets, observed)
        write_output(records, args.output)
        any_yes = any(present == "yes" for _, present in records)
        sys.exit(0 if any_yes else 2)

    except Exception as exc:
        logging.exception("Internal error: %s", exc)
        sys.exit(1)


if __name__ == "__main__":
    main()
