#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Usage:
    validate_taxa --taxa <RANK_TAXA.csv> --taxdump <NEWTAXDUMP>
    validate_taxa --taxa <RANK_TAXA.csv> --taxdump <NEWTAXDUMP> --gtdb-map \
                    --gtdb-info /path/to/dir/with_gtdb_info \
                    [--out mapping.csv] [--threshold 20.0]

Purpose:
    1) Validate that each (rank, taxa) pair in a taxa CSV exists in the given
       NCBI taxdump (taxdump.json), and that its rank matches.

       Additionally, if the taxa string is GTDB-prefixed (one of:
         d__, p__, c__, o__, f__, g__, s__
       ), treat it as a GTDB taxon (not an NCBI taxon) and validate it against:
         <gtdb-info>/gtdb_r226.dic

    2) (Optional: --gtdb-map) Resolve each NCBI input to its NCBI phylum
       (using ONLY taxdump.json parents) and map that NCBI phylum
       to GTDB phylum(s) using the Phylum sheets in:
        - ncbi_vs_gtdb_bacteria.xlsx
        - ncbi_vs_gtdb_archaea.xlsx

       For GTDB-prefixed inputs, mapping columns are left blank:
         ncbi_phylum == "" and gtdb_phylum == ""

Inputs:
    -t, --taxa      : CSV with header 'rank,taxa'
    -d, --taxdump   : Directory containing 'taxdump.json'
    --gtdb-info     : Directory containing GTDB resources:
                      - 'gtdb_r226.dic'
                      - 'ncbi_vs_gtdb_bacteria.xlsx'
                      - 'ncbi_vs_gtdb_archaea.xlsx'

Notes:
    - (rank,taxa) pairs are deduplicated
    - Rank normalisation: 'realm' and 'domain' are treated as 'superkingdom'
    - GTDB mapping is only for Bacteria/Archaea; others get blank gtdb_phylum

Author: Akito Shima (asuq)
Email: asuq.4096@gmail.com
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import re
import sys
from pathlib import Path
from typing import Any

import pandas as pd


logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

if sys.version_info < (3, 12):
    logging.fatal(
        f"Python 3.12 or newer is required. Current version: {sys.version.split()[0]}"
    )
    sys.exit(1)


TAXDUMP_TIMESTAMP: str = "20240914"
GTDB_RELEASE: str = "R226"

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
TAXA_PCT_RE = re.compile(r"([^\s%]+)\s*([0-9]+(?:\.[0-9]+)?)\s*%?")

VIRUSES_REALM = {
    "adnaviria",
    "duplodnaviria",
    "monodnaviria",
    "riboviria",
    "ribozyviria",
    "singelaviria",
    "varidnaviria",
}

GTDB_PREFIX_TO_NORM_RANK: dict[str, str] = {
    "d__": "superkingdom",
    "p__": "phylum",
    "c__": "class",
    "o__": "order",
    "f__": "family",
    "g__": "genus",
    "s__": "species",
}
GTDB_PREFIXES: tuple[str, ...] = tuple(GTDB_PREFIX_TO_NORM_RANK.keys())


def is_gtdb_taxon(taxa: str) -> bool:
    tx = taxa.strip()
    return any(tx.startswith(pfx) for pfx in GTDB_PREFIXES)


def gtdb_expected_norm_rank(taxa: str) -> str | None:
    tx = taxa.strip()
    for pfx, rk in GTDB_PREFIX_TO_NORM_RANK.items():
        if tx.startswith(pfx):
            return rk
    return None


def load_gtdb_dic(dic_path: Path) -> set[str]:
    """
    Load gtdb_r226.dic (ASCII, one taxon per line) into a set.

    Entries are kept case-sensitive and include the GTDB prefix (e.g. 'p__Bacillota').
    Blank lines and lines starting with '#' are ignored.
    """
    taxa: set[str] = set()
    try:
        with dic_path.open("r", encoding="ascii") as fh:
            for line in fh:
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                taxa.add(s)
    except UnicodeDecodeError as e:
        raise ValueError(f"{dic_path} must be ASCII text: {e}") from e

    return taxa


def _fmt_not_found_ncbi(taxa: str, rank: str) -> str:
    return (
        f"Taxon '{taxa}' not found in the NCBI Taxonomy ({TAXDUMP_TIMESTAMP}) "
        f"(requested rank={rank})."
    )


def _fmt_not_found_gtdb(taxa: str, rank: str, dic_path: Path) -> str:
    return (
        f"GTDB taxon '{taxa}' not found in {dic_path.name} ({GTDB_RELEASE}) "
        f"(requested rank={rank})."
    )


def _fmt_rank_mismatch_ncbi(
    taxa: str, requested_norm_rank: str, found_ranks: list[str]
) -> str:
    found = ", ".join(found_ranks) if found_ranks else "unknown"
    return (
        f"Taxon '{taxa}' found but rank mismatch: taxdump has {found}; "
        f"requested '{requested_norm_rank}'. "
        "Note: this tool normalises 'realm'/'domain' to 'superkingdom'. "
        "Use the canonical rank from the pinned snapshot."
    )


def _fmt_rank_mismatch_gtdb(
    taxa: str, requested_norm_rank: str, expected_norm_rank: str
) -> str:
    return (
        f"GTDB taxon '{taxa}' found but rank mismatch: prefix implies '{expected_norm_rank}', "
        f"requested '{requested_norm_rank}'. "
        "Note: this tool normalises 'realm'/'domain' to 'superkingdom'."
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate taxa CSV against taxdump.json. Optionally map to GTDB phyla.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
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
        "--gtdb-map",
        action="store_true",
        help="Also map NCBI inputs to GTDB phylum using the Phylum sheet (bacteria/archaea).",
    )
    parser.add_argument(
        "--gtdb-info",
        "--ncbi_to_gtdb",
        dest="gtdb_info",
        type=Path,
        help=(
            "Directory containing GTDB resources: "
            "'gtdb_r226.dic', 'ncbi_vs_gtdb_bacteria.xlsx', 'ncbi_vs_gtdb_archaea.xlsx'. "
        ),
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("mapping.csv"),
        help="Output CSV for mapping results (header: rank,taxa,ncbi_phylum,gtdb_phylum).",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=20.0,
        help="Percent threshold (strictly greater) for reporting GTDB phyla.",
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
        raise FileNotFoundError(f"File not found: {path}")
    if kind == "dir" and not path.is_dir():
        raise FileNotFoundError(f"Directory not found: {path}")
    if valid_exts:
        if path.suffix.lower() not in tuple(e.lower() for e in valid_exts):
            raise ValueError(
                f"Invalid extension for {path.name}. "
                f"Expected one of: {', '.join(valid_exts)}"
            )


# --------------------------------------------------------------------------- #
#  Taxdump loader (JSON-only)
# --------------------------------------------------------------------------- #
class TaxdumpIndex:
    """
    Case-insensitive name -> list[(taxid, normalised_rank)] built from taxdump.json.

    Expected top-level keys:
      - names     : {taxid: "canonical name"}
      - ranks     : {taxid: "rank"}  (e.g., "domain", "phylum", "genus"...)
      - ancestors : {taxid: {rank: int}}
                    (e.g. ancestors["1775750"]["phylum"] == 3055124)

    Rank normalisation:
      - 'domain' → 'superkingdom'
      - 'realm'  → 'superkingdom'
    """

    _RANK_NORMALISE = {"domain": "superkingdom", "realm": "superkingdom"}

    def __init__(self, taxdump_json: Path) -> None:
        self._name_index: dict[str, list[tuple[str, str]]] = {}
        self._ranks_by_id: dict[str, str] = {}
        self._names_by_id: dict[str, str] = {}
        self.ancestors: dict[str, dict[str, int]] = {}
        self._load(taxdump_json)

    @staticmethod
    def _norm_rank(rank: str | None) -> str:
        if not rank:
            return ""
        r = rank.strip().lower()
        return TaxdumpIndex._RANK_NORMALISE.get(r, r)

    def _load(self, path: Path) -> None:
        with path.open("r", encoding="utf-8") as fh:
            data: Any = json.load(fh)

        if not isinstance(data, dict) or "names" not in data or "ranks" not in data:
            raise ValueError("taxdump.json must be a dict with keys 'names' and 'ranks'")

        names: dict[str, str] = data["names"]  # {taxid: name}
        ranks: dict[str, str] = data["ranks"]  # {taxid: rank}
        ancestors_in: dict[str, dict[str, int]] = data.get("ancestors", {})

        for tid, name in names.items():
            if not name:
                continue
            rank = self._norm_rank(ranks.get(tid, ""))
            key = name.strip().lower()
            if key == "viruses" and rank != "superkingdom":
                # Treat “Viruses” as superkingdom for legacy compatibility
                rank = "superkingdom"
            self._name_index.setdefault(key, []).append((str(tid), rank))
            self._ranks_by_id[str(tid)] = rank
            self._names_by_id[str(tid)] = name.strip()

        anc: dict[str, dict[str, int]] = {}
        for tid, table in ancestors_in.items():
            table_norm: dict[str, int] = {}
            for rk, val in table.items():
                # val should be int (0, positive taxid, or negative sentinel)
                table_norm[self._norm_rank(rk)] = (
                    int(val) if isinstance(val, (int, float)) else 0
                )
            anc[str(tid)] = table_norm
        self.ancestors = anc

    def lookup(self, name: str) -> list[tuple[str, str]]:
        return list(self._name_index.get(name.strip().lower(), []))

    def resolve_to_phylum(
        self, name: str, requested_rank: str
    ) -> tuple[str | None, str | None]:
        """
        Resolve 'name' (at given requested_rank) to:
          (ncbi_phylum_name, superkingdom_name) via ancestors table.
        Returns (None, None) if the phylum cannot be determined.
        """
        entries = self.lookup(name)
        if not entries:
            return None, None

        req_rank = INPUT_RANK_TO_COLUMN.get(requested_rank, requested_rank)
        candidates = [tid for (tid, rk) in entries if rk == req_rank]
        if not candidates:
            # Fallback to any entry if the requested rank doesn't match any
            candidates = [tid for (tid, _) in entries]
            logging.warning(
                f"Name {name} found but no entry matched requested rank {requested_rank}; using first match"
            )

        # Prefer canonical name match; else take the first
        chosen: str = candidates[0]
        if len(candidates) > 1:
            name_lc = name.strip().lower()
            for tid in candidates:
                canon = self._names_by_id.get(tid, "")
                if canon and canon.strip().lower() == name_lc:
                    chosen = tid
                    break
            else:
                logging.warning(f"Ambiguous matches for {name}; picked taxid {chosen}")

        # Ancestor table lookup
        anc = self.ancestors.get(chosen, {})

        # Determine phylum taxid:
        ph_val = anc.get("phylum", 0)
        ph_tid: str | None = None
        if ph_val > 0:
            ph_tid = str(ph_val)
        elif ph_val == 0 and self._ranks_by_id.get(chosen) == "phylum":
            # if this node is itself a phylum
            ph_tid = chosen
        else:
            # negative value means the node is above the phylum level (no phylum to assign)
            ph_tid = None

        # Determine superkingdom taxid:
        sk_val = anc.get("superkingdom", 0)
        sk_tid: str | None = None
        if sk_val > 0:
            sk_tid = str(sk_val)
        elif sk_val == 0 and self._ranks_by_id.get(chosen) == "superkingdom":
            sk_tid = chosen

        ph_name = self._names_by_id.get(ph_tid) if ph_tid else None
        sk_name = self._names_by_id.get(sk_tid) if sk_tid else None
        return ph_name, sk_name


# --------------------------------------------------------------------------- #
#  Core logic
# --------------------------------------------------------------------------- #
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

    seen: set[tuple[str, str]] = set()
    dedup: list[tuple[str, str]] = []
    for p in pairs:
        if p not in seen:
            seen.add(p)
            dedup.append(p)
    return dedup


# --------------------------------------------------------------------------- #
#  GTDB crosswalk helpers (Phylum sheet only)
# --------------------------------------------------------------------------- #
def parse_gtdb_list_field(field: Any) -> list[tuple[str, float]]:
    """
    Parse strings like:
      "p__Pseudomonadota 99.51%, p__Bacteroidota 0.08%, ..."
    Return: [(GTDB_phylum_plain, percent_float), ...]
    """
    out: list[tuple[str, float]] = []
    if field is None:
        return out

    s = str(field)
    if not s or s.strip().lower() in {"nan", "none"}:
        return out

    parts = [p.strip() for p in s.split(",") if str(p).strip()]
    for part in parts:
        m = TAXA_PCT_RE.search(part)
        if not m:
            continue
        name_raw, pct_raw = m.group(1), m.group(2)
        # strip p__ prefix if present; output should be plain names
        name = name_raw[3:] if name_raw.startswith("p__") else name_raw

        # safe float parse: keep only digits and a single dot
        digits: list[str] = []
        saw_dot = False
        for ch in pct_raw:
            if ch == "." and not saw_dot:
                digits.append(ch)
                saw_dot = True
            elif ch.isdigit():
                digits.append(ch)
        pct = float("".join(digits)) if digits else 0.0
        out.append((name.strip(), pct))
    return out


def load_phylum_sheet(xlsx_path: Path) -> dict[str, list[tuple[str, float]]]:
    """
    Load only the 'Phylum' sheet and return mapping:
      ncbi_phylum_plain_lower -> [(gtdb_phylum_plain, percent_float), ...]
    """
    validate_file(xlsx_path, "file", (".xlsx",))
    df = pd.read_excel(xlsx_path, sheet_name="Phylum")
    cols = {str(c).strip().lower(): c for c in df.columns}
    if "ncbi phylum" not in cols or "list of r226 phyla" not in cols:
        raise ValueError(
            f"Phylum sheet in {xlsx_path} missing columns 'NCBI phylum' / 'List of R226 phyla'"
        )

    mapping: dict[str, list[tuple[str, float]]] = {}
    for _, row in df.iterrows():
        phylum_val = row[cols["ncbi phylum"]]
        if phylum_val is None:
            continue
        ncbi_phylum = str(phylum_val).strip()
        if not ncbi_phylum or ncbi_phylum.lower() in {"nan", "none"}:
            continue
        if ncbi_phylum.startswith("p__"):
            ncbi_phylum = ncbi_phylum[3:]
        list_field = row[cols["list of r226 phyla"]]
        parsed = parse_gtdb_list_field(list_field)
        mapping[ncbi_phylum.strip().lower()] = parsed
    return mapping


# --------------------------------------------------------------------------- #
#  Main
# --------------------------------------------------------------------------- #
def main() -> int:
    args = parse_args()

    # Validate inputs & extensions
    validate_file(args.taxa, "file", (".csv",))
    validate_file(args.taxdump, "dir")
    taxdump_json = args.taxdump / "taxdump.json"
    validate_file(taxdump_json, "file", (".json",))

    selections = load_taxa_file(args.taxa)
    logging.info(f"{len(selections)} selection(s) loaded from --taxa.")

    any_gtdb = any(is_gtdb_taxon(taxa) for _, taxa in selections)
    if (args.gtdb_map or any_gtdb) and not args.gtdb_info:
        logging.error(
            "GTDB resources are required (input contains GTDB-prefixed taxa and/or --gtdb-map was set). "
            "Provide --gtdb-info pointing to a directory containing gtdb_r226.dic "
            "(and the ncbi_vs_gtdb_*.xlsx files if using --gtdb-map)."
        )
        return 1

    gtdb_dic_path: Path | None = None
    gtdb_taxa_set: set[str] = set()
    if args.gtdb_info:
        validate_file(args.gtdb_info, "dir")
        gtdb_dic_path = args.gtdb_info / "gtdb_r226.dic"
        if any_gtdb:
            validate_file(gtdb_dic_path, "file")
            logging.info(f"Loading {gtdb_dic_path} for GTDB taxon validation")
            gtdb_taxa_set = load_gtdb_dic(gtdb_dic_path)

    logging.info(f"Loading taxdump.json from {taxdump_json}")
    taxidx = TaxdumpIndex(taxdump_json)

    errors: list[str] = []
    for rank, taxa in selections:
        requested_norm_rank = INPUT_RANK_TO_COLUMN[rank]

        if is_gtdb_taxon(taxa):
            if not gtdb_dic_path:
                errors.append(
                    f"GTDB taxon '{taxa}' present but --gtdb-info was not provided."
                )
                continue

            if taxa not in gtdb_taxa_set:
                errors.append(_fmt_not_found_gtdb(taxa, rank, gtdb_dic_path))
                continue

            expected_norm = gtdb_expected_norm_rank(taxa) or ""
            if expected_norm and expected_norm != requested_norm_rank:
                errors.append(
                    _fmt_rank_mismatch_gtdb(taxa, requested_norm_rank, expected_norm)
                )
            continue

        entries = taxidx.lookup(taxa)
        if not entries:
            errors.append(_fmt_not_found_ncbi(taxa, rank))
            continue

        ranks_found = sorted({rk for _, rk in entries if rk})
        if not any(rk == requested_norm_rank for rk in ranks_found):
            errors.append(_fmt_rank_mismatch_ncbi(taxa, requested_norm_rank, ranks_found))

    if errors:
        logging.error("Validation failed:")
        for msg in errors:
            logging.error(f"  {msg}")
        logging.error(
            f"NCBI names must match the taxonomy snapshot ({TAXDUMP_TIMESTAMP}) used to build the GTDB ncbi_vs_gtdb tables."
        )
        return 1

    logging.info("All taxa validated successfully.")

    if not args.gtdb_map:
        return 0

    if not args.gtdb_info:
        logging.error(
            "--gtdb-map requires --gtdb-info pointing to a directory with the GTDB resources."
        )
        return 1

    bact_xlsx = args.gtdb_info / "ncbi_vs_gtdb_bacteria.xlsx"
    arch_xlsx = args.gtdb_info / "ncbi_vs_gtdb_archaea.xlsx"

    # Load mappings if present; warn if absent
    bact_map: dict[str, list[tuple[str, float]]] = {}
    arch_map: dict[str, list[tuple[str, float]]] = {}
    if bact_xlsx.is_file():
        bact_map = load_phylum_sheet(bact_xlsx)
    else:
        logging.warning(f"Missing {bact_xlsx}; bacterial mappings will be empty.")
    if arch_xlsx.is_file():
        arch_map = load_phylum_sheet(arch_xlsx)
    else:
        logging.warning(f"Missing {arch_xlsx}; archaeal mappings will be empty.")

    # Prepare output
    out_path = args.out
    with out_path.open("w", newline="", encoding="utf-8") as out_fh:
        writer = csv.writer(out_fh, delimiter=",", lineterminator="\n")
        writer.writerow(["rank", "taxa", "ncbi_phylum", "gtdb_phylum"])

        threshold = float(args.threshold)

        for rank, taxa in selections:
            if is_gtdb_taxon(taxa):
                # For GTDB inputs, leave mapping columns blank.
                writer.writerow([rank, taxa, "", ""])
                continue

            ncbi_phylum, superkingdom = taxidx.resolve_to_phylum(taxa, rank)

            # 2) Skip Eukaryotes / Viruses (but keep row).
            sk_lower = (superkingdom or "").strip().lower()
            if sk_lower in {"eukaryota", "eukarya"}:
                logging.info(f"{taxa} is Eukaryotes")
                writer.writerow([rank, taxa, ncbi_phylum or "", ""])
                continue
            if sk_lower == "viruses" or sk_lower in VIRUSES_REALM:
                logging.info(f"{taxa} is Viruses")
                writer.writerow([rank, taxa, ncbi_phylum or "", ""])
                continue

            # If cannot resolve phylum, still emit row (empty gtdb_phylum)
            if ncbi_phylum is None:
                writer.writerow([rank, taxa, "", ""])
                logging.info(f"no ncbi_phylum found for {taxa}")
                continue

            # 3) Choose the correct crosswalk
            phylum_key = ncbi_phylum.strip()
            if phylum_key.startswith("p__"):
                phylum_key = phylum_key[3:]
            key_lower = phylum_key.lower()

            mapping: dict[str, list[tuple[str, float]]] | None = None
            if sk_lower.startswith("bac") or sk_lower == "bacteria":
                mapping = bact_map
            elif sk_lower.startswith("arc") or sk_lower == "archaea":
                mapping = arch_map
            else:
                # fallback: pick the dictionary that contains the phylum if any
                if key_lower in bact_map and key_lower not in arch_map:
                    mapping = bact_map
                elif key_lower in arch_map and key_lower not in bact_map:
                    mapping = arch_map
                elif key_lower in bact_map and key_lower in arch_map:
                    mapping = bact_map  # deterministic fallback
                else:
                    mapping = None

            if not mapping:
                logging.warning(
                    f"No GTDB mapping source available for NCBI phylum '{ncbi_phylum}' "
                    f"(taxon '{taxa}', superkingdom='{superkingdom or ''}'); emitting empty gtdb_phylum."
                )
                writer.writerow([rank, taxa, ncbi_phylum, ""])
                continue

            candidate = mapping.get(key_lower)
            if not candidate:
                logging.warning(
                    f"NCBI phylum '{ncbi_phylum}' (taxon '{taxa}') not present in selected GTDB mapping; "
                    f"emitting empty gtdb_phylum."
                )
                writer.writerow([rank, taxa, ncbi_phylum, ""])
                continue

            # 4) Filter strictly > threshold; emit one row per passing GTDB phylum
            passed = [g for (g, pct) in candidate if pct > threshold]
            if not passed:
                logging.warning(
                    f"No GTDB phylum above threshold {threshold}% for NCBI phylum '{ncbi_phylum}' "
                    f"(taxon '{taxa}'); emitting empty gtdb_phylum."
                )
                writer.writerow([rank, taxa, ncbi_phylum, ""])
                continue

            for gtdb_phylum in passed:
                writer.writerow([rank, taxa, ncbi_phylum, gtdb_phylum])

    logging.info(f"Mapping written to {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
