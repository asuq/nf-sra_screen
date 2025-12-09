#!/usr/bin/env python3
"""
Check presence of target GTDB taxa (or GTDB phylum names) in a single-sample taxonomy TSV.

Exit status:
0 : at least one target is present (i.e., some row has present == 'yes')
1 : internal error (I/O problems, malformed input, unexpected exceptions)
2 : no target is present (i.e., all rows have present == 'no' or no rows)

Usage:
    python check_singlem_taxa.py \
        --input sample_taxonomy.tsv \
        --gtdb_taxa-list gtdb_taxa_to_check.txt \
        --output gtdb_taxa_presence.tsv \
        --verbose

Input TSV:
    -i, --input : TSV output from singlem pipe (header includes 'taxonomy')
      Example taxonomy cell:
        'Root; d__Bacteria; p__Bacillota; c__Bacilli; ...'

Targets list:
    -g, --gtdb_taxa-list : Text file with a header line 'gtdb_taxa',
      followed by one target per line.

      Each target line can be either:
        (A) GTDB-prefixed taxon: d__/p__/c__/o__/f__/g__/s__  (e.g. 'g__Bacillus')
            -> matched against *any* token in the taxonomy path at any rank.
        (B) Unprefixed GTDB phylum name (e.g. 'Bacillota')
            -> treated as a phylum target and matched against the taxonomy phylum token.

      Matching is case-insensitive. For GTDB-prefixed targets, the prefix is preserved
      and included in matching (i.e., rank-aware).

Output:
    -o, --output : TSV with columns (gtdb_taxa, present)
      where 'present' is 'yes' or 'no' for each requested target label.

Required Packages: pandas

Author: Akito Shima (asuq)
Email: asuq.4096@gmail.com
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
from collections.abc import Iterable, Sequence
from pathlib import Path

import pandas as pd


if sys.version_info < (3, 12):
    print(
        f"Python 3.12 or newer is required. Current version: %s",
        sys.version.split()[0],
    )
    sys.exit(1)


GTDB_PREFIXES: tuple[str, ...] = ("d__", "p__", "c__", "o__", "f__", "g__", "s__")
GTDB_PREFIX_RE = re.compile(r"^(d__|p__|c__|o__|f__|g__|s__)", flags=re.IGNORECASE)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse CLI arguments."""
    p = argparse.ArgumentParser(
        description=(
            "Check whether a single-sample taxonomy TSV contains target taxa. "
            "Targets can be GTDB-prefixed taxa (d__/p__/c__/o__/f__/g__/s__) or "
            "unprefixed GTDB phylum names (treated as phylum targets)."
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
        "-g",
        "--gtdb_taxa-list",
        type=Path,
        required=True,
        help="Path to targets list text file (header 'gtdb_taxa', one per line).",
    )
    p.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Path to output TSV (columns: gtdb_taxa, present).",
    )
    p.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose (DEBUG) logging.",
    )
    return p.parse_args(argv)


def setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(levelname)s: %(message)s")


def normalise_target_label(label: str) -> str | None:
    """
    Normalise a target label into a canonical key for matching.

    Rules:
      1) If the label begins with a GTDB prefix (d__/p__/c__/o__/f__/g__/s__),
         keep the prefix and normalise to: '<prefix><name>' all lower-case.
         Example: 'G__Bacillus' -> 'g__bacillus'
      2) Otherwise, treat the label as an unprefixed GTDB phylum name and
         normalise to: 'p__<name>' all lower-case.
         Example: 'Bacillota' -> 'p__bacillota'

    Returns None for empty/whitespace-only labels.
    """
    s = label.strip()
    if not s:
        return None

    m = GTDB_PREFIX_RE.match(s)
    if m:
        prefix = m.group(1).lower()
        rest = s[len(m.group(1)) :].strip()
        if not rest:
            return None
        return f"{prefix}{rest.lower()}"

    # Unprefixed labels are interpreted as phylum names (legacy behaviour).
    return f"p__{s.lower()}"


def extract_normalised_tokens_from_taxonomy(taxonomy: str) -> set[str]:
    """
    Extract all GTDB-prefixed tokens from a taxonomy path and normalise them.

    Example input:
        'Root; d__Bacteria; p__Bacillota; c__Bacilli'
    Output set contains:
        {'d__bacteria', 'p__bacillota', 'c__bacilli'}

    Tokens without a recognised prefix are ignored.
    """
    if not isinstance(taxonomy, str):
        return set()

    # Primary delimiter is ';' for SingleM; be permissive if files use '|'.
    raw = taxonomy
    parts = raw.split(";") if ";" in raw else raw.split("|")

    out: set[str] = set()
    for token in (t.strip() for t in parts):
        if not token:
            continue
        m = GTDB_PREFIX_RE.match(token)
        if not m:
            continue
        prefix = m.group(1).lower()
        rest = token[len(m.group(1)) :].strip()
        if not rest:
            continue
        out.add(f"{prefix}{rest.lower()}")
    return out


def read_observed_taxa_from_input(input_path: Path) -> set[str]:
    """
    Read the input TSV and return a set of observed GTDB tokens (normalised keys).

    Only the 'taxonomy' column is used.
    """
    logging.info("Reading taxonomy TSV: %s", input_path)
    df = pd.read_csv(input_path, sep="\t", dtype=str, usecols=["taxonomy"])

    observed: set[str] = set()
    for tax in df["taxonomy"].dropna():
        observed |= extract_normalised_tokens_from_taxonomy(tax)

    logging.info("Observed %d unique GTDB tokens in file.", len(observed))
    logging.debug("Observed tokens (normalised): %s", sorted(observed))
    return observed


def read_target_labels(list_path: Path) -> list[tuple[str, str]]:
    """
    Read targets list file and return pairs of (original_label, normalised_key).

    Header:
      - First non-empty line must be 'gtdb_taxa'.
    """
    logging.info("Reading targets list: %s", list_path)
    raw_text = list_path.read_text(encoding="utf-8")

    targets: list[tuple[str, str]] = []
    for i, raw in enumerate(raw_text.splitlines()):
        line = raw.strip().lstrip("\ufeff")
        if not line:
            continue

        if i == 0 and line.lower() == "gtdb_taxa":
            continue

        key = normalise_target_label(line)
        if key is None:
            continue
        targets.append((line, key))

    if not targets:
        logging.warning("No targets found in list file after processing header.")
    else:
        logging.info("Loaded %d target(s).", len(targets))
        logging.debug("Targets (original -> normalised): %s", targets)
    return targets


def compute_presence_records(
    targets: Sequence[tuple[str, str]], observed: set[str]
) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    for original, key in targets:
        present = "yes" if key in observed else "no"
        records.append((original, present))
    return records


def write_output(records: Iterable[tuple[str, str]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df_out = pd.DataFrame(records, columns=["gtdb_taxa", "present"])
    df_out.to_csv(output_path, sep="\t", index=False)
    logging.info("Wrote results: %s", output_path)


def main(argv: Sequence[str] | None = None) -> None:
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
        if not args.gtdb_taxa_list.exists():
            raise FileNotFoundError(f"Targets list not found: {args.gtdb_taxa_list}")

        observed = read_observed_taxa_from_input(args.input)
        targets = read_target_labels(args.gtdb_taxa_list)
        records = compute_presence_records(targets, observed)
        write_output(records, args.output)

        any_yes = any(present == "yes" for _, present in records)
        sys.exit(0 if any_yes else 2)

    except Exception as exc:
        logging.exception("Internal error: %s", exc)
        sys.exit(1)


if __name__ == "__main__":
    main()
