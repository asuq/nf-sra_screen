#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Combined utilities for working with NCBI taxdump snapshots used by GTDB.

Purpose:
  1) make-taxidlineage: Create new_taxdump-style taxidlineage.dmp from an old-format taxdump.
  2) jsonify-taxdump (default action): Build taxdump.json using nodes/names/taxidlineage.

Behaviour:
  - Default CLI (positional TAXDUMP dir): If taxidlineage.dmp(.gz) is absent,
    first generate it from nodes.dmp(.gz) with exclude-self & skip-deleted,
    then write taxdump.json (or taxdump.json.gz if lineage is gz).
    If lineage exists, directly write taxdump.json (compression matches lineage).

Attribution:
  Portions of the JSON taxonomy logic are adapted from BlobToolKit (Sanger Tree of Life; sanger-tol).
  See: https://github.com/blobtoolkit

Required Packages: ujson

Author: Akito Shima (asuq)
Email: asuq.4096@gmail.com
"""

import argparse
import csv
import gzip
import io
import logging
import re
import sys
from collections.abc import Iterable
from pathlib import Path
from typing import TextIO, cast

import ujson


logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

if sys.version_info < (3, 12):
    logging.fatal(
        f"Python 3.12 or newer is required. Current version: {sys.version.split()[0]}"
    )
    sys.exit(1)

# =============================================================================
#  PART 1: make_taxidlineage
# =============================================================================
# ---- Type aliases ----
type TaxID = int
type ParentMap = dict[TaxID, TaxID]
type ChildrenMap = dict[TaxID, list[TaxID]]

ALIAS_TO_SUPERKINGDOM = {
    "realm": "superkingdom",
    "domain": "superkingdom",
    "superkingdom": "superkingdom",
}

# -------------------------- generic helpers --------------------------


def open_text(path: Path, mode: str = "rt", encoding: str = "utf-8") -> TextIO:
    """
    Open a plain or gzipped text file, based on suffix.
    Supports 'rt' and 'wt' modes.
    """
    if "t" not in mode:
        raise ValueError("open_text only supports text modes ('rt' or 'wt').")
    if ".gz" in path.suffixes:
        return cast(TextIO, gzip.open(path, mode=mode, encoding=encoding))
    return cast(TextIO, path.open(mode=mode, encoding=encoding))


def parse_dmp_line(line: str) -> list[str]:
    """
    Parse NCBI .dmp line with <tab>|<tab> separators and trailing <tab>|.
    Return fields with surrounding whitespace stripped and trailing empty removed.
    """
    parts: list[str] = [p.strip() for p in line.split("|")]
    while parts and parts[-1] == "":
        parts.pop()
    return parts


def write_ncbi_dmp_line(outfh: TextIO, taxid: TaxID, lineage_tokens: Iterable[TaxID]) -> None:
    """
    Write: tax_id <TAB>|<TAB> <space-separated lineage> <space if non-empty> <TAB>|

    Example (non-empty lineage):
      2798934\t|\t131567 2157 ... 2798928 \t|
    Example (empty lineage):
      1\t|\t\t|
    """
    tokens_list: list[str] = [str(x) for x in lineage_tokens]
    lineage_str: str = (" ".join(tokens_list) + " ") if tokens_list else ""
    outfh.write(f"{taxid}\t|\t{lineage_str}\t|\n")


def find_file(directory: Path, base: str) -> Path | None:
    """
    Return directory/base if it exists else directory/(base + '.gz') if it exists, else None.
    """
    p = directory / base
    if p.exists():
        return p
    gz = directory / f"{base}.gz"
    if gz.exists():
        return gz
    return None


# -------------------------- parsing --------------------------


def read_nodes(nodes_path: Path) -> tuple[ParentMap, ChildrenMap]:
    """
    Parse nodes.dmp(.gz) into:
      parent_of[taxid] = parent_taxid
      children_of[parent_taxid] = sorted list of child taxids (no self-loops)
    """
    parent_of: ParentMap = {}
    children_of: ChildrenMap = {}

    with open_text(nodes_path, "rt") as fh:
        for i, raw in enumerate(fh, 1):
            line: str = raw.strip()
            if not line:
                continue
            cols: list[str] = parse_dmp_line(line)
            if len(cols) < 2:
                raise ValueError(f"Malformed nodes.dmp line {i}: {line!r}")
            try:
                tax_id: TaxID = int(cols[0])
                parent_id: TaxID = int(cols[1])
            except ValueError as e:
                raise ValueError(f"Non-integer at nodes.dmp line {i}: {cols[:2]}") from e

            parent_of[tax_id] = parent_id

            # Skip self-parent edges (e.g., 1 -> 1) to avoid repeating root lines.
            if tax_id != parent_id:
                children_of.setdefault(parent_id, []).append(tax_id)

    # Deterministic traversal order
    for _, kids in children_of.items():
        kids.sort()

    return parent_of, children_of


def read_merged(merged_path: Path | None) -> dict[TaxID, TaxID]:
    """
    Parse merged.dmp(.gz) into a map old_taxid -> new_taxid.
    """
    merged: dict[TaxID, TaxID] = {}
    if merged_path is None or not merged_path.exists():
        return merged

    with open_text(merged_path, "rt") as fh:
        for i, raw in enumerate(fh, 1):
            line: str = raw.strip()
            if not line:
                continue
            cols: list[str] = parse_dmp_line(line)
            if len(cols) < 2:
                raise ValueError(f"Malformed merged.dmp line {i}: {line!r}")
            try:
                old_id: TaxID = int(cols[0])
                new_id: TaxID = int(cols[1])
            except ValueError as e:
                raise ValueError(f"Non-integer at merged.dmp line {i}: {cols[:2]}") from e
            merged[old_id] = new_id
    return merged


def read_delnodes(delnodes_path: Path | None) -> set[TaxID]:
    """
    Parse delnodes.dmp(.gz) into a set of deleted/retired taxids.
    """
    deleted: set[TaxID] = set()
    if delnodes_path is None or not delnodes_path.exists():
        return deleted

    with open_text(delnodes_path, "rt") as fh:
        for i, raw in enumerate(fh, 1):
            line: str = raw.strip()
            if not line:
                continue
            cols: list[str] = parse_dmp_line(line)
            if not cols:
                continue
            try:
                deleted.add(int(cols[0]))
            except ValueError as e:
                raise ValueError(f"Non-integer at delnodes.dmp line {i}: {cols[:1]}") from e
    return deleted


def find_roots(parent_of: ParentMap) -> list[TaxID]:
    """
    Roots are nodes where parent_taxid == taxid or where the parent is missing.
    Typically the only true root is 1, but we handle orphans defensively.
    """
    roots: list[TaxID] = []
    for t, p in parent_of.items():
        if p == t or p not in parent_of:
            roots.append(t)
    roots.sort()
    return roots


# -------------------------- lineage writing --------------------------


def dfs_write_taxidlineage(
    parent_of: ParentMap,
    children_of: ChildrenMap,
    outfh: TextIO,
    include_root_id: bool = False,
    include_self: bool = True,
) -> None:
    """
    Stream taxidlineage rows by iterative DFS from each root.
    """
    roots: list[TaxID] = find_roots(parent_of)
    visited: set[TaxID] = set()

    for root in roots:
        stack: list[tuple[TaxID, int]] = [(root, 0)]  # (node, next_child_index)
        path: list[TaxID] = [root]  # root -> ... -> node
        path_set: set[TaxID] = {root}  # for O(1) cycle checks

        def current_lineage_tokens() -> list[TaxID]:
            tokens: list[TaxID] = list(path)  # root..node
            if not include_root_id and tokens and tokens[0] == 1:
                tokens = tokens[1:]
            if not include_self and tokens:
                tokens = tokens[:-1]
            return tokens

        if root not in visited:
            write_ncbi_dmp_line(outfh, root, current_lineage_tokens())
            visited.add(root)

        while stack:
            node, idx = stack[-1]
            kids: list[TaxID] = children_of.get(node, [])
            if idx < len(kids):
                child: TaxID = kids[idx]
                stack[-1] = (node, idx + 1)  # advance index

                # Cycle guard: skip if child already on current path
                if child in path_set:
                    continue

                stack.append((child, 0))
                path.append(child)
                path_set.add(child)

                if child not in visited:
                    write_ncbi_dmp_line(outfh, child, current_lineage_tokens())
                    visited.add(child)
            else:
                stack.pop()
                popped: TaxID = path.pop()
                path_set.discard(popped)

    # Fallback for any stray nodes not reached from declared roots.
    if len(visited) < len(parent_of):
        for t in sorted(set(parent_of) - visited):
            chain: list[TaxID] = []
            cur: TaxID = t
            seen: set[TaxID] = set()
            while cur in parent_of and cur not in seen:
                chain.append(cur)
                seen.add(cur)
                p: TaxID = parent_of[cur]
                if p == cur:
                    break
                cur = p
            chain.reverse()
            if not include_root_id and chain and chain[0] == 1:
                chain = chain[1:]
            if not include_self and chain:
                chain = chain[:-1]
            write_ncbi_dmp_line(outfh, t, chain)


def build_taxidlineage(
    nodes_path: Path,
    out_path: Path,
    include_root_id: bool = False,
    include_self: bool = True,
    merged_path: Path | None = None,
    add_merged: bool = False,
    delnodes_path: Path | None = None,
    skip_deleted: bool = False,
) -> None:
    """
    Orchestrate reading inputs and writing the output file. All files may be .gz.
    """
    parent_of, children_of = read_nodes(nodes_path)

    with open_text(out_path, "wt") as outfh:
        dfs_write_taxidlineage(
            parent_of=parent_of,
            children_of=children_of,
            outfh=outfh,
            include_root_id=include_root_id,
            include_self=include_self,
        )

        if add_merged:
            merged = read_merged(merged_path)
            deleted = read_delnodes(delnodes_path) if skip_deleted else set()

            lineage_cache: dict[TaxID, list[TaxID]] = {}

            def lineage_of(taxid: TaxID) -> list[TaxID]:
                if taxid in lineage_cache:
                    return lineage_cache[taxid]
                chain: list[TaxID] = []
                cur: TaxID = taxid
                seen: set[TaxID] = set()
                while cur in parent_of and cur not in seen:
                    chain.append(cur)
                    seen.add(cur)
                    p: TaxID = parent_of[cur]
                    if p == cur:
                        break
                    cur = p
                chain.reverse()
                if not include_root_id and chain and chain[0] == 1:
                    chain = chain[1:]
                if not include_self and chain:
                    chain = chain[:-1]
                lineage_cache[taxid] = chain
                return chain

            for old_id in sorted(merged):
                if skip_deleted and old_id in deleted:
                    continue
                new_id: TaxID = merged[old_id]
                if new_id not in parent_of:
                    # Replacement not present in this snapshot; skip.
                    continue
                write_ncbi_dmp_line(outfh, old_id, lineage_of(new_id))


def main_taxidlineage(argv: list[str] | None = None) -> int:
    """
    CLI entry for the lineage writer.
    """
    parser = argparse.ArgumentParser(
        prog="make-taxidlineage",
        description="Create taxidlineage.dmp (new_taxdump-style) from old-format taxdump (nodes.dmp).",
    )
    parser.add_argument(
        "--taxdump",
        required=True,
        type=Path,
        help="Directory containing old-format taxdump (must include nodes.dmp or nodes.dmp.gz).",
    )
    parser.add_argument(
        "--out",
        "-o",
        type=Path,
        default=Path("taxidlineage.dmp"),
        help="Output path; '.gz' => gzip. Default: taxidlineage.dmp",
    )
    parser.add_argument(
        "--exclude-self",
        action="store_true",
        default=False,
        help="Exclude the node’s own TaxID from the lineage. (default: include self)",
    )
    parser.add_argument(
        "--add-merged",
        action="store_true",
        default=False,
        help="Also add rows for obsolete taxids using merged.dmp(.gz).",
    )
    parser.add_argument(
        "--skip-deleted",
        action="store_true",
        default=False,
        help="When --add-merged, skip obsolete IDs listed in delnodes.dmp(.gz).",
    )

    args = parser.parse_args(argv)
    taxdump_dir = Path(args.taxdump)
    if args.out.exists():
        logging.warning(f"{args.out} already exists; skipping")
        return 0

    # Locate nodes and optional merged/delnodes (gz-aware)
    nodes_path = find_file(taxdump_dir, "nodes.dmp")
    if nodes_path is None:
        logging.error(f"nodes.dmp(.gz) not found in {taxdump_dir}")
        return 2

    merged_path: Path | None = find_file(taxdump_dir, "merged.dmp") if args.add_merged else None
    delnodes_path: Path | None = (
        find_file(taxdump_dir, "delnodes.dmp") if args.skip_deleted else None
    )
    if args.add_merged and merged_path is None:
        logging.warning(
            f"--add-merged set but merged.dmp(.gz) not found; merged IDs will be skipped."
        )

    try:
        logging.info(f"Building {args.out} from nodes.dmp(.gz)")
        build_taxidlineage(
            nodes_path=nodes_path,
            out_path=args.out,
            include_root_id=False,
            include_self=not args.exclude_self,
            merged_path=merged_path,
            add_merged=args.add_merged and (merged_path is not None),
            delnodes_path=delnodes_path,
            skip_deleted=args.skip_deleted,
        )
    except Exception as e:
        logging.error(f"Failed to build lineage: {e}")
        return 2

    return 0


# =============================================================================
#  PART 2: jsonify_taxdump
#  Attribution: adapted from BlobToolKit (Sanger Tree of Life; sanger-tol).
# =============================================================================.


def stream_file(filename):
    """
    Stream a file, line by line.

    Automatically detect gzipped files based on suffix.
    """
    path = filename if isinstance(filename, Path) else Path(str(filename))
    if ".gz" in path.suffixes:
        try:
            return gzip.open(path, "rt")
        except OSError:
            return None
    try:
        return path.open("r", encoding="utf-8")
    except OSError:
        return None


def write_file(filename, data, plain=False):
    """
    Write a file, use suffix to determine type and compression.

    - types: '.json'
    - compression: None, '.gz'

    write_file('variable.json.gz')
    """
    fname = str(filename)
    if ".json" in fname:
        content = ujson.dumps(data, indent=1, escape_forward_slashes=False)
    elif filename == "STDOUT":
        sys.stdout.write(ujson.dumps(data, indent=1, escape_forward_slashes=False) + "\n")
        return True
    elif filename == "STDERR":
        sys.stderr.write(ujson.dumps(data, indent=1, escape_forward_slashes=False) + "\n")
        return True
    elif plain:
        content = "\n".join(data)
    elif ".csv" in fname or ".tsv" in fname:
        output = io.StringIO()
        if ".csv" in fname:
            writer = csv.writer(output, quoting=csv.QUOTE_NONNUMERIC)
        else:
            writer = csv.writer(output, delimiter="\t")
        for row in data:
            writer.writerow(row)
        content = output.getvalue()
    else:
        content = data
    if ".gz" in fname:
        try:
            with gzip.open(fname, "wt") as fh:
                fh.write(content)
        except OSError:
            return False
    else:
        try:
            Path(fname).write_text(content, encoding="utf-8")
        except OSError:
            return False
    return True


class Taxdump:
    """Class for working with NCBI taxonomy."""

    __slots__ = ["directory", "ancestors", "names", "ranks"]

    def __init__(self, directory, **kwargs):
        """Init Taxdump class."""
        self.directory = directory if isinstance(directory, Path) else Path(directory)
        self.ancestors = {}
        self.ranks = {}
        self.names = {}
        if kwargs:
            self.update_data(**kwargs)
        else:
            self.load_ranks()
            self.load_names()
            self.load_ancestors()

    def update_data(self, **kwargs):
        """Update values and keys for an existing field."""
        for key, value in kwargs.items():
            setattr(self, key, {int(taxid): data for taxid, data in value.items()})

    def _file(self, base: str) -> Path:
        """
        Locate base or base.gz under self.directory.
        """
        p = self.directory / base
        if p.exists():
            return p
        gz = self.directory / f"{base}.gz"
        if gz.exists():
            return gz
        raise FileNotFoundError(f"Required file not found: {base}(.gz) in {self.directory}")

    def load_ranks(self):
        """Load ranks from file."""
        try:
            filename = self._file("nodes.dmp")
            for line in stream_file(filename):
                row = self.parse_taxdump_row(line)
                if len(row) > 2:
                    self.ranks[int(row[0])] = row[2]
        except TypeError:
            logging.error(f"Unable to parse {filename}.")
            exit(1)

    def load_names(self):
        """Load names from file."""
        try:
            filename = self._file("names.dmp")
            for line in stream_file(filename):
                row = self.parse_taxdump_row(line)
                if len(row) > 3 and row[3] == "scientific name":
                    self.names[int(row[0])] = row[1]
        except TypeError:
            logging.error(f"Unable to parse {filename}.")
            exit(1)

    def load_ancestors(self):
        """Load ancestors from file."""
        try:
            filename = self._file("taxidlineage.dmp")
            saw_domain = False
            saw_realm = False

            ranks_list = self.list_ranks()

            for line in stream_file(filename):
                row = self.parse_taxdump_row(line)
                if len(row) > 1 and row[1]:
                    taxid = int(row[0])

                    mapping: dict[str, int] = {}
                    for tok in row[1].split():
                        try:
                            anc_id = int(tok)
                        except ValueError:
                            continue
                        if anc_id not in self.ranks:
                            continue

                        raw = self.ranks[anc_id]
                        if raw == "domain":
                            saw_domain = True
                        elif raw == "realm":
                            saw_realm = True

                        key = ALIAS_TO_SUPERKINGDOM.get(raw, raw)
                        # First-wins policy for each key; do not overwrite if already set.
                        if key in ranks_list and key not in mapping:
                            mapping[key] = anc_id

                    # Consider the node's own rank as well (after ancestors).
                    if taxid in self.ranks:
                        raw_self = self.ranks[taxid]
                        if raw_self == "domain":
                            saw_domain = True
                        elif raw_self == "realm":
                            saw_realm = True

                        key_self = ALIAS_TO_SUPERKINGDOM.get(raw_self, raw_self)
                        if key_self in ranks_list and key_self not in mapping:
                            mapping[key_self] = taxid

                    # Fill fixed rank axis with negative sentinels, top-down.
                    last = 0
                    for rank in ranks_list:
                        if rank in mapping:
                            last = -mapping[rank]
                        else:
                            mapping[rank] = last

                    self.ancestors[taxid] = mapping

            # Warn once if we observed modern top ranks being aliased.
            if saw_domain or saw_realm:
                kinds = []
                if saw_domain:
                    kinds.append("domain")
                if saw_realm:
                    kinds.append("realm")
                logging.warning(
                    f"Detected {kinds} in taxonomy; mapping to 'superkingdom' in JSON 'ancestors' for compatibility."
                )

        except TypeError:
            logging.error(f"Unable to parse {filename}.")
            exit(1)

    def values_to_dict(self):
        """Create a dict of values."""
        data = {}
        for key in ("ancestors", "names", "ranks"):
            if hasattr(self, key):
                data[key] = getattr(self, key)
        return data

    def lineage(self, taxid):
        """Create a dict of rank-name pairs for a taxid."""
        lineages = {}
        try:
            ancestors = self.ancestors[taxid]
        except (KeyError, ValueError):
            return {}
        for rank in self.list_ranks():
            if rank in ancestors and ancestors[rank] > 0:
                lineages.update({rank: self.names[ancestors[rank]]})
        return lineages

    @staticmethod
    def parse_taxdump_row(line):
        """Parse an ncbi taxdump file."""
        return re.split(r"\s*\|\s*", line)[:-1]

    @staticmethod
    def list_ranks():
        """Return a list of taxonomic ranks."""
        return [
            "superkingdom",
            "kingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species",
        ]


# =============================================================================
#  Pipeline (default action)
# =============================================================================


def build_lineage_and_json(taxdump_dir: Path) -> int:
    """
    If taxidlineage.dmp(.gz) is missing, generate it with exclude-self & skip-deleted
    from nodes.dmp(.gz), then write taxdump.json (or taxdump.json.gz if lineage .gz).
    If lineage exists, directly write JSON.
    """
    logging.info(f"Checking if taxidlineage.dmp exists")

    taxdump_dir = Path(taxdump_dir)

    # Prefer uncompressed lineage, then gz
    lineage_path = find_file(taxdump_dir, "taxidlineage.dmp")
    if lineage_path is None:
        # Build lineage: locate inputs (gz-aware)
        logging.info(f"taxidlineage.dmp doesn't exist, building from nodes.dmp")
        nodes = find_file(taxdump_dir, "nodes.dmp")
        if nodes is None:
            logging.error(f"nodes.dmp(.gz) not found in {taxdump_dir}")
            return 2

        merged = find_file(taxdump_dir, "merged.dmp")  # optional
        delnodes = find_file(taxdump_dir, "delnodes.dmp")  # optional

        # Generate uncompressed lineage by default
        lineage_path = taxdump_dir / "taxidlineage.dmp"
        try:
            build_taxidlineage(
                nodes_path=nodes,
                out_path=lineage_path,
                include_root_id=False,  # exclude root '1'
                include_self=False,  # --exclude-self
                merged_path=merged,
                add_merged=False,  # spec did not request --add-merged
                delnodes_path=delnodes,
                skip_deleted=True,  # --skip-deleted
            )
        except Exception as e:  # surface helpful message
            logging.error(f"Failed to build lineage: {e}")
            return 2

    json_out = taxdump_dir / (
        "taxdump.json.gz" if ".gz" in lineage_path.suffixes else "taxdump.json"
    )
    if json_out.exists():
        logging.warning(f"{json_out} already exists; skipping")
        return 0

    # Build JSON using Taxdump class (gz-aware readers)
    try:
        logging.info(f"Building taxdump.json from .dmp files")
        taxdump = Taxdump(taxdump_dir)
        ok = write_file(json_out, taxdump.values_to_dict())
        if not ok:
            logging.error(f"Failed to write {json_out}")
            return 2
    except Exception as e:
        logging.error(f"Failed to write JSON: {e}")
        return 2

    logging.info(f"Job completed")
    return 0


# =============================================================================
#  Main
# =============================================================================


def main() -> int:
    # If the first argument looks like the lineage tool, dispatch there; otherwise
    # run the default pipeline (build lineage if needed, then write JSON file).
    if len(sys.argv) > 1 and sys.argv[1] in {
        "make-taxidlineage",
        "make_taxidlineage",
        "taxidlineage",
    }:
        return main_taxidlineage(sys.argv[2:])
    if len(sys.argv) < 2:
        sys.stderr.write(
            "Usage: jsonify_taxdump.py <taxdump_dir>  OR  jsonify_taxdump.py make-taxidlineage [options]\n"
        )
        return 2
    taxdump_dir = Path(sys.argv[1])
    return build_lineage_and_json(taxdump_dir)


if __name__ == "__main__":
    sys.exit(main())
