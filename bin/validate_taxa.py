#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Usage:
    validate_taxa --taxa <RANK_TAXA.csv> --taxdump <NEWTAXDUMP_DIR>

Purpose:
    Validate that each (rank, taxa) pair in a taxa CSV exists in the given
    BlobToolKit-style new_taxdump (taxdump.json), and that its rank matches.
    Raises an error (exit 1) if any entry is invalid.

Inputs:
    -t, --taxa      : CSV with header 'rank,taxa'
    -d, --taxdump   : Directory containing 'taxdump.json' (BlobToolKit style)

Notes:
    - (rank,taxa) pairs are deduplicated (set semantics)
    - Rank normalization: 'realm' and 'domain' are treated as 'superkingdom'

Author: Akito Shima (ASUQ)
Email: asuq.4096@gmail.com
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import sys
from pathlib import Path
from typing import Any


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
    p = argparse.ArgumentParser(
        description="Validate taxa CSV against BlobToolKit taxdump.json.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--taxa", "-t", type=Path, required=True, help="CSV with header 'rank,taxa'."
    )
    p.add_argument(
        "--taxdump",
        "-d",
        type=Path,
        required=True,
        help="Directory containing 'taxdump.json'.",
    )
    return p.parse_args()


def validate_file(
    path: Path, kind: str, valid_exts: tuple[str, ...] | None = None
) -> None:
    if kind == "file" and not path.is_file():
        raise FileNotFoundError(f"File not found: {path}")
    if kind == "dir" and not path.is_dir():
        raise FileNotFoundError(f"Directory not found: {path}")
    if valid_exts and path.suffix.lower() not in tuple(e.lower() for e in valid_exts):
        raise ValueError(
            f"Invalid extension for {path.name}. Expected one of: {', '.join(valid_exts)}"
        )


class TaxdumpIndex(object):
    """
    Case-insensitive name → list[(taxid, normalized_rank)] built from BlobToolKit-style taxdump.json.

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
            data: Any = json.load(fh)

        if not isinstance(data, dict) or "names" not in data or "ranks" not in data:
            raise ValueError(
                "taxdump.json must be a dict with keys 'names' and 'ranks' as per BlobToolKit"
            )

        names = data["names"]  # {taxid: name}
        ranks = data["ranks"]  # {taxid: rank}

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
        return list(self._map.get(name.strip().lower(), []))


def load_taxa_file(taxa_csv: Path) -> list[tuple[str, str]]:
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

    return list(set(pairs))


def main() -> int:
    args = parse_args()
    validate_file(args.taxa, "file", (".csv",))
    validate_file(args.taxdump, "file", (".json",))

    logging.info(f"Loading taxdump.json from {args.taxdump}")
    taxidx = TaxdumpIndex(args.taxdump)

    selections = load_taxa_file(args.taxa)
    logging.info(f"{len(selections)} selection(s) loaded from --taxa.")

    errors: list[str] = []
    for rank_req, taxa_name in selections:
        # realm/domain → superkingdom
        requested_norm_rank = INPUT_RANK_TO_COLUMN[rank_req]
        entries = taxidx.lookup(taxa_name)

        if not entries:
            errors.append(f"Taxon '{taxa_name}' not found (requested rank={rank_req})")
            continue

        ranks_found = sorted({rk for _, rk in entries if rk})
        if not any(rk == requested_norm_rank for rk in ranks_found):
            errors.append(
                f"Taxon '{taxa_name}' found but rank mismatch: found {ranks_found or 'unknown'}, expected '{requested_norm_rank}' (from requested '{rank_req}')"
            )

    if errors:
        logging.error("Validation failed:")
        for msg in errors:
            logging.error("  " + msg)
        sys.exit(1)

    logging.info("All taxa validated successfully.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
