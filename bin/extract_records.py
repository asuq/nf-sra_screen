#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Usage:
    extract_records --blobtable <BLOBTABLE.csv> --fasta <ASSEMBLY.fasta> \
                    --taxa <RANK_TAXA.csv> --taxdump <NEWTAXDUMP> \
                    [--out_dir <OUT_DIR>]

Purpose:
    Extract and segregate FASTA records per (rank, taxa) using BlobTools output
    and new_taxdump (taxdump.json).

Inputs:
    -b, --blobtable : CSV output from blobtools2 with columns:
        identifiers,gc,length,ncount,coverage,superkingdom,kingdom,phylum,class,order,family,genus,species
    -f, --fasta     : Assembly FASTA (record.id must match 'identifiers')
    -t, --taxa      : CSV with header 'rank,taxa' (one or more rows)
    -d, --taxdump   : Directory containing 'taxdump.json'. No .dmp parsing.
    -o, --out_dir   : Output directory path [default: ./]

Notes:
    - (rank,taxa) pairs are deduplicated
    - If a taxon name is not found in taxdump.json → WARN (continue)
    - If a taxon exists but at a different rank → WARN (continue)
    - Rank normalization: 'realm' and 'domain' are treated as 'superkingdom'
    - If a selection yields 0 identifiers, no files are written

Required Packages: Biopython, Pandas

Author: Akito Shima (ASUQ)
Email: asuq.4096@gmail.com
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import pandas as pd
from Bio import SeqIO


logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


if sys.version_info < (3, 10):
    logging.fatal(
        f"Python 3.10 or newer is required. Current version: {sys.version.split()[0]}"
    )
    sys.exit(1)


INPUT_RANK_TO_COLUMN: dict[str, str] = {
    "realm": "superkingdom",
    "domain": "superkingdom",
    "superkingdom": "superkingdom",
    "kingdom": "kingdom",
    "phylum": "phylum",
    "class": "class",
    "order": "order",
    "family": "family",
    "genus": "genus",
    "species": "species",
}

ALLOWED_RANKS: tuple[str, ...] = tuple(INPUT_RANK_TO_COLUMN.keys())


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract FASTA records per (rank,taxa) using BlobTools CSV and taxdump.json.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--blobtable",
        "-b",
        type=Path,
        required=True,
        help="BlobTools table CSV (must contain 'identifiers' and standard rank columns).",
    )
    parser.add_argument(
        "--fasta",
        "-f",
        type=Path,
        required=True,
        help="Assembly FASTA whose headers match 'identifiers'.",
    )
    parser.add_argument(
        "--taxa",
        "-t",
        type=Path,
        required=True,
        help="CSV with header 'rank,taxa' listing taxa to extract.",
    )
    parser.add_argument(
        "--taxdump",
        "-d",
        type=Path,
        required=True,
        help="Directory containing 'taxdump.json'.",
    )
    parser.add_argument(
        "--out_dir",
        "-o",
        type=Path,
        default=Path("./"),
        help="Output directory for flat files.",
    )
    return parser.parse_args()


def validate_file(
    path: Path, kind: str = "file", valid_exts: tuple[str, ...] | None = None
) -> None:
    """
    Validate that the path exists, matches the expected kind, and
    (optionally) has one of the allowed extensions.
    """

    if kind == "file" and not path.is_file():
        raise FileNotFoundError(f"Required file not found: {path}")
    if kind == "dir" and not path.is_dir():
        raise FileNotFoundError(f"Required directory not found: {path}")

    if valid_exts:
        if path.suffix.lower() not in tuple(e.lower() for e in valid_exts):
            raise ValueError(
                f"Invalid extension for '{path.name}'. "
                f"Expected one of: {', '.join(valid_exts)}"
            )


def sanitize_filename_component(s: str) -> str:
    """Replace non-alphanumeric/._- chars with underscores."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", s.strip())


# --------------------------------------------------------------------------- #
#  Taxdump loader (JSON-only)
# --------------------------------------------------------------------------- #
class TaxdumpIndex(object):
    """
    Case-insensitive name → list[(taxid, normalized_rank)] built from taxdump.json.

    Expected top-level keys:
      - names : {taxid: "canonical name"}
      - ranks : {taxid: "rank"}  (e.g., "domain", "phylum", "genus"...)

    Rank normalization:
      - 'domain' → 'superkingdom'
      - 'realm'  → 'superkingdom'
    """

    _RANK_NORMALIZE = {"domain": "superkingdom", "realm": "superkingdom"}

    def __init__(self, taxdump_json: Path) -> None:
        self._map: dict[str, list[tuple[str, str]]] = {}
        self._load(taxdump_json)

    @staticmethod
    def _norm_rank(rank: str | None) -> str:
        if not rank:
            return ""
        r = rank.strip().lower()
        return TaxdumpIndex._RANK_NORMALIZE.get(r, r)

    def _load(self, path: Path) -> None:
        with path.open("r", encoding="utf-8") as fh:
            data = json.load(fh)

        if not isinstance(data, dict) or "names" not in data or "ranks" not in data:
            raise ValueError(
                "taxdump.json must be a dict with keys 'names' and 'ranks'"
            )

        names = data["names"]  # {taxid: name}
        ranks = data["ranks"]  # {taxid: rank}

        # Build name → [(taxid, norm_rank)]
        index: dict[str, list[tuple[str, str]]] = {}
        for tid, name in names.items():
            if not name:
                continue
            rank_raw = ranks.get(tid, "")
            rank = self._norm_rank(rank_raw)
            key = str(name).strip().lower()
            index.setdefault(key, []).append((str(tid), rank))

        self._map = index

    def lookup(self, name: str) -> list[tuple[str, str]]:
        """Return list of (taxid, normalized_rank) for a given (case-insensitive) name."""
        return list(self._map.get(name.strip().lower(), []))


# --------------------------------------------------------------------------- #
#  Core logic
# --------------------------------------------------------------------------- #
def load_taxa_file(taxa_csv: Path) -> list[tuple[str, str]]:
    """Read --taxa CSV (header rank,taxa), deduplicate preserving order."""
    with taxa_csv.open("r", newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None or {"rank", "taxa"} - set(reader.fieldnames):
            raise ValueError("--taxa must have header 'rank,taxa'")

        pairs: list[tuple[str, str]] = []
        for row in reader:
            rank = (row["rank"] or "").strip().lower()
            taxa = (row["taxa"] or "").strip()
            if not rank or not taxa:
                continue
            if rank not in ALLOWED_RANKS:
                raise ValueError(
                    f"Invalid rank '{rank}'. Allowed: {', '.join(ALLOWED_RANKS)}"
                )
            pairs.append((rank, taxa))

    seen, dedup = set(), []
    for p in pairs:
        if p not in seen:
            seen.add(p)
            dedup.append(p)
    return dedup


def select_identifiers(df: pd.DataFrame, rank: str, taxa: str) -> set[str]:
    col = INPUT_RANK_TO_COLUMN[rank]  # normalize 'realm'/'domain' → 'superkingdom'
    nulls = {"", "na", "none", "no-hit", "-", "nan"}
    series = df[col].astype(str).str.strip().str.lower()
    mask = (~series.isin(nulls)) & (series.eq(taxa.strip().lower()))
    return set(df.loc[mask, "identifiers"].astype(str))


def dispatch_fasta(fasta: Path, dispatch: dict[str, list[Path]]) -> None:
    """
    Stream FASTA once, writing matching records to relevant output files.

    dispatch: identifiers -> list of Path for output FASTA files
    """
    handles: dict[Path, Any] = {}
    try:
        for rec in SeqIO.parse(str(fasta), "fasta"):
            targets = dispatch.get(rec.id)
            if not targets:
                continue
            for target in targets:
                if target not in handles:
                    handles[target] = target.open("w")
                SeqIO.write(rec, handles[target], "fasta")
    finally:
        for h in handles.values():
            h.close()


# --------------------------------------------------------------------------- #
#  Main
# --------------------------------------------------------------------------- #
def main() -> int:
    args = parse_args()

    # Validate inputs & extensions
    validate_file(args.blobtable, "file", (".csv",))
    validate_file(args.taxa, "file", (".csv",))
    validate_file(args.fasta, "file", (".fasta", ".fa", ".fna", ".fas", ".mpfa"))
    validate_file(args.taxdump, "dir")
    taxdump_json = args.taxdump / "taxdump.json"
    validate_file(taxdump_json, "file", (".json",))

    logging.info(f"Loading taxdump.json from {taxdump_json}")
    taxidx = TaxdumpIndex(taxdump_json)

    selections = load_taxa_file(args.taxa)
    logging.info(f"{len(selections)} selection(s) loaded from --taxa.")

    # Validate taxa presence/rank in taxdump (warn-only)
    for rank, taxa in selections:
        entries = taxidx.lookup(taxa)
        if not entries:
            logging.warning(f"Taxon '{taxa}' not found in taxdump.json (rank={rank})")
            continue

        expected = INPUT_RANK_TO_COLUMN[rank].lower()  # realm/domain → superkingdom
        ranks_found = sorted({(rk or "").lower() for _, rk in entries if rk})
        if expected not in ranks_found:
            logging.warning(
                f"Taxon '{taxa}' found but ranks differ ({sorted(ranks_found) or 'unknown'}), expected '{expected}'"
            )

    # Load blobtable and validate columns
    df = pd.read_csv(args.blobtable, dtype={"identifiers": "string"})
    need_cols = {"identifiers", *{INPUT_RANK_TO_COLUMN[r] for r in ALLOWED_RANKS}}
    missing = need_cols - set(df.columns)
    if missing:
        sys.exit(f"--blobtable missing columns: {', '.join(sorted(missing))}")

    args.out_dir.mkdir(parents=True, exist_ok=True)

    dispatch: dict[str, list[Path]] = defaultdict(list)
    summary_rows: list[tuple[str, str, int, str, str]] = []

    # For each (rank,taxa), write IDs CSV and prep FASTA target
    for rank, taxa in selections:
        ids = select_identifiers(df, rank, taxa)
        if not len(ids):
            logging.warning(f"{rank}={taxa}: no hits")
            summary_rows.append((rank, taxa, 0, "", ""))
            continue

        fname = sanitize_filename_component(f"{rank}_{taxa}")
        ids_path = args.out_dir / f"{fname}.ids.csv"
        with ids_path.open("w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["identifiers"])
            for ident in sorted(ids):
                w.writerow([ident])
        logging.info(f"{rank}={taxa}: {len(ids)} id(s) -> {ids_path.name}")

        fa_path = args.out_dir / f"{fname}.fasta"
        for ident in ids:
            dispatch[ident].append(fa_path)

        summary_rows.append((rank, taxa, len(ids), str(ids_path), str(fa_path)))

    # Stream FASTA once and emit records to their respective outputs
    if dispatch:
        logging.info(f"Streaming FASTA: {args.fasta.name}")
        dispatch_fasta(args.fasta, dispatch)
    else:
        logging.warning("No hits across all selections")

    # Write overview summary.csv
    summary_csv = args.out_dir / "summary.csv"
    with summary_csv.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["rank", "taxa", "n_contigs", "output_ids_csv", "output_fasta"])
        w.writerows(summary_rows)
    logging.info(f"Wrote summary: {summary_csv}")

    logging.info("All done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
