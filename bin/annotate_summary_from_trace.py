#!/usr/bin/env python3
"""
Annotate summary.tsv with scheduler-level failure reasons derived from trace.tsv.

Usage
-----
    annotate_summary_from_trace.py SUMMARY_TSV TRACE_TSV

Arguments
---------
SUMMARY_TSV
    The summary.tsv produced by run_append_summary.sh

TRACE_TSV
    Nextflow trace file

Behaviour
---------
- Only the last attempt of each Nextflow task is considered, using
  the largest numeric task_id for each `name`.
- Existing rows in summary.tsv are annotated in-place (note column).
- If a (sra, srr) has a failing last-attempt task in trace.tsv but no
  row in summary.tsv, a new summary row is created using metadata from:
    <outdir>/metadata/<sra>/<sra>.filtered.csv
    <outdir>/metadata/<sra>/<sra>.skipped.csv

The metadata CSVs are expected to have at least these headers:
    accession,run_accession,instrument_platform,instrument_model,
    library_source,library_strategy,assembler
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable

from tqdm import tqdm


LOG = logging.getLogger("annotate_summary")


# --------------------------------------------------------------------------- #
# CLI and logging
# --------------------------------------------------------------------------- #
def setup_logging(verbosity: int) -> None:
    """Configure root logger based on requested verbosity level."""
    level = logging.WARNING
    if verbosity == 1:
        level = logging.INFO
    elif verbosity >= 2:
        level = logging.DEBUG

    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Annotate summary.tsv with scheduler-level failure reasons "
            "derived from a Nextflow trace file."
        ),
    )
    parser.add_argument(
        "summary_tsv",
        type=Path,
        help="Path to summary.tsv produced by run_append_summary.sh",
    )
    parser.add_argument(
        "trace_txt",
        type=Path,
        help="Path to Nextflow trace file (tab-separated)",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Increase log verbosity (can be repeated)",
    )
    return parser.parse_args(argv)


# --------------------------------------------------------------------------- #
# Core classification / parsing helpers
# --------------------------------------------------------------------------- #
def classify_failure(
    process_name: str,
    status: str | None,
    exit_code: int | None,
    native_id: str | None,
) -> str | None:
    """
    Map (status, exit_code, native_id) with error message

    Categories:
    - job submission error
    - time limit error
    - memory / resource limit error
    - generic execution error
    """
    status_norm = (status or "").strip().upper()
    native_norm = (native_id or "").strip()

    if status_norm in {"", "COMPLETED", "CACHED"}:
        return None

    if status_norm == "CANCELLED":
        return f"{process_name}: cancelled (pipeline aborted or upstream failure)"

    if status_norm == "ABORTED":
        return f"{process_name}: aborted by Nextflow (upstream error)"

    # Job submission error: FAILED but no real scheduler job id
    if status_norm == "FAILED" and (
        not native_norm or native_norm in {"-", "0", "NA", "N/A", "NULL", "null"}
    ):
        if exit_code is not None:
            return f"{process_name}: job submission error " f"(no native id; exit {exit_code})"
        return f"{process_name}: job submission error (no native id; unknown exit code)"

    # From here on, assume the job actually ran on the executor.
    if exit_code is None:
        return f"{process_name}: execution error (status {status_norm}, unknown exit)"

    # Memory / resource limits – killed by signal.
    # 137 = 128 + 9   (SIGKILL) -> often OOM
    # 139 = SIGSEGV   -> crash / memory bug
    # 140 = 128 + 12  (SIGUSR2; cluster-dependent)
    if exit_code in {137, 139, 140}:
        return (
            f"{process_name}: memory / resource limit error "
            f"(killed by signal; exit {exit_code})"
        )

    # Time limit – many schedulers use SIGTERM for walltime.
    # 143 = 128 + 15  (SIGTERM)
    if exit_code == 143:
        return f"{process_name}: time limit error (terminated by SIGTERM; exit 143)"

    # Environment / container issues.
    if exit_code == 127:
        return (
            f"{process_name}: execution error "
            "(command not found; exit 127 – check container/environment)"
        )

    # Generic non-zero exit.
    return f"{process_name}: execution error (exit {exit_code})"


def parse_tag_from_name(name: str | None) -> tuple[str | None, str | None]:
    """
    Extract (sra, srr) from a task name string.

    Expected format:
        'PROCESS_NAME (SRAX:SRR12345)'

    Returns
    -------
    (sra, srr) or (None, None) if parsing fails.
    """
    if not name:
        return None, None

    stripped = name.strip()
    if "(" not in stripped or not stripped.endswith(")"):
        return None, None

    tag = stripped[stripped.rfind("(") + 1 : -1]  # content inside parentheses
    if ":" not in tag:
        return None, None

    sra, srr = tag.split(":", 1)
    sra = sra.strip()
    srr = srr.strip()
    return (sra or None), (srr or None)


# --------------------------------------------------------------------------- #
# Summary loading / writing
# --------------------------------------------------------------------------- #
def load_summary(summary_path: Path) -> tuple[list[dict[str, str]], list[str]]:
    """
    Load summary.tsv and ensure it has a 'note' column.

    Returns
    -------
    rows, fieldnames
    """
    LOG.info("Loading summary.tsv from %s", summary_path)

    with summary_path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows: list[dict[str, str]] = list(reader)
        fieldnames: list[str] = list(reader.fieldnames or [])

    if "note" not in fieldnames:
        LOG.debug("Adding missing 'note' column to summary.tsv")
        fieldnames.append("note")
        for row in rows:
            row["note"] = ""

    return rows, fieldnames


def build_index_by_key(
    rows: Iterable[dict[str, str]],
) -> dict[tuple[str, str], int]:
    """
    Build an index mapping (sra, srr) -> row index.

    Missing values are treated as empty strings.
    """
    index: dict[tuple[str, str], int] = {}
    for i, row in enumerate(rows):
        sra = row.get("sra", "") or ""
        srr = row.get("srr", "") or ""
        index[(sra, srr)] = i
    return index


def write_summary(
    summary_path: Path,
    rows: list[dict[str, str]],
    fieldnames: list[str],
) -> None:
    """Write updated summary.tsv atomically via a temporary file."""
    tmp_path = summary_path.with_suffix(summary_path.suffix + ".tmp")

    LOG.info(
        "Writing updated summary to temporary file %s (final: %s)",
        tmp_path,
        summary_path,
    )

    with tmp_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    tmp_path.replace(summary_path)
    LOG.info("Successfully replaced %s with updated summary.tsv", summary_path)


# --------------------------------------------------------------------------- #
# Trace handling (last attempt per task)
# --------------------------------------------------------------------------- #
def iter_trace_records(trace_path: Path) -> Iterable[dict[str, str]]:
    """Yield trace records as dicts from a tab-separated file."""
    with trace_path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        yield from reader


def last_attempt_records(trace_path: Path) -> list[dict[str, str]]:
    """
    Select the last attempt per Nextflow task using the largest task_id.

    Grouping key is the full `name` field, which is typically:
        'PROCESS_NAME (TAG)'
    where TAG encodes SRA/SRR in your pipeline.
    """
    LOG.info("Selecting last attempt per task using task_id")
    best_by_name: dict[str, tuple[int, dict[str, str]]] = {}

    for record in tqdm(
        iter_trace_records(trace_path),
        desc="Scanning trace",
        unit="task",
    ):
        name = (record.get("name") or "").strip()
        if not name:
            continue

        task_id_str = (record.get("task_id") or "").strip()
        try:
            task_id = int(task_id_str)
        except ValueError:
            LOG.debug("Skipping record with non-integer task_id=%r", task_id_str)
            continue

        prev = best_by_name.get(name)
        if prev is None or task_id > prev[0]:
            best_by_name[name] = (task_id, record)

    LOG.info("Found %d distinct tasks (by name)", len(best_by_name))
    return [rec for _, rec in best_by_name.values()]


def collect_reasons_from_trace(
    trace_path: Path,
) -> dict[tuple[str, str], list[str]]:
    """
    Parse trace.tsv and return a mapping (sra, srr) -> list of messages.

    Only the last attempt (largest task_id) of each task is considered.
    """
    last_records = last_attempt_records(trace_path)
    reasons: dict[tuple[str, str], list[str]] = defaultdict(list)

    for record in tqdm(
        last_records,
        desc="Collecting failure reasons",
        unit="task",
    ):
        status = (record.get("status") or "").strip()
        if not status:
            continue

        status_upper = status.upper()
        if status_upper in {"COMPLETED", "CACHED"}:
            continue

        name = (record.get("name") or "").strip()
        native_id = record.get("native_id") or ""
        exit_str = (record.get("exit") or "").strip()

        try:
            exit_code: int | None = int(exit_str) if exit_str else None
        except ValueError:
            exit_code = None

        sra, srr = parse_tag_from_name(name)
        if not sra or not srr:
            # e.g. VALIDATE_TAXA or untagged processes
            continue

        process_name = name.split()[0] if name else "UNKNOWN"
        message = classify_failure(process_name, status, exit_code, native_id)

        if message:
            reasons[(sra, srr)].append(message)

    LOG.info(
        "Collected scheduler reasons for %d (sra, srr) combinations",
        len(reasons),
    )
    return reasons


# --------------------------------------------------------------------------- #
# Metadata lookup from ./output/metadata
# --------------------------------------------------------------------------- #
def load_metadata_for_sras(
    metadata_root: Path,
    sras: set[str],
) -> dict[tuple[str, str], dict[str, str]]:
    """
    Build an index of metadata rows keyed by (accession, run_accession).

    It looks for, per SRA:

      metadata_root/<sra>/<sra>.filtered.csv
      metadata_root/<sra>/<sra>.skipped.csv

    Returns
    -------
    mapping: (sra, srr) -> metadata_row_dict
    """
    mapping: dict[tuple[str, str], dict[str, str]] = {}

    if not sras:
        return mapping

    LOG.info(
        "Loading metadata for %d SRA(s) from %s",
        len(sras),
        metadata_root,
    )

    for sra in sorted(sras):
        candidates = [
            metadata_root / sra / f"{sra}.filtered.csv",
            metadata_root / sra / f"{sra}.skipped.csv",
        ]

        for csv_path in candidates:
            if not csv_path.is_file():
                continue

            LOG.debug("Reading metadata from %s", csv_path)

            with csv_path.open(newline="") as fh:
                reader = csv.DictReader(fh)
                for row in reader:
                    acc = (row.get("accession") or "").strip()
                    run = (row.get("run_accession") or "").strip()
                    if not run:
                        continue

                    key_sra = acc or sra
                    key = (key_sra, run)
                    # Last one wins if duplicates exist; should be identical anyway.
                    mapping[key] = row

    LOG.info("Built metadata index with %d (sra, srr) keys", len(mapping))
    return mapping


# --------------------------------------------------------------------------- #
# Annotation logic
# --------------------------------------------------------------------------- #
def annotate_summary(
    summary_rows: list[dict[str, str]],
    fieldnames: list[str],
    reasons: dict[tuple[str, str], list[str]],
    metadata_root: Path,
) -> None:
    """
    Merge scheduler reasons into the 'note' column per (sra, srr).

    If a (sra, srr) is not present in summary_rows, a new row is created
    using metadata from metadata_root (if available).
    """
    if not reasons:
        LOG.info("No scheduler failures to annotate; summary.tsv unchanged")
        return

    index = build_index_by_key(summary_rows)
    LOG.debug("Built index for %d summary rows", len(index))

    # Figure out which SRAs we need metadata for (only those missing from summary).
    missing_keys = [key for key in reasons if key not in index]
    missing_sras: set[str] = {sra for sra, _ in missing_keys}

    metadata_index: dict[tuple[str, str], dict[str, str]] = {}
    if missing_sras:
        metadata_index = load_metadata_for_sras(metadata_root, missing_sras)
    else:
        LOG.info("No missing (sra, srr) in summary; metadata lookup not needed")

    for (sra, srr), messages in reasons.items():
        unique_messages = sorted(set(messages))
        scheduler_note = "; ".join(unique_messages)

        row_idx = index.get((sra, srr))
        if row_idx is None:
            LOG.info(
                "No existing summary row for (%s, %s); creating a new one",
                sra,
                srr,
            )
            meta_row = metadata_index.get((sra, srr), {})

            new_row: dict[str, str] = {}
            for field in fieldnames:
                if field == "sra":
                    new_row[field] = sra
                elif field == "srr":
                    new_row[field] = srr
                elif field == "platform":
                    new_row[field] = (meta_row.get("instrument_platform") or "").strip()
                elif field == "model":
                    new_row[field] = (meta_row.get("instrument_model") or "").strip()
                elif field == "strategy":
                    new_row[field] = (meta_row.get("library_strategy") or "").strip()
                elif field == "assembler":
                    new_row[field] = (meta_row.get("assembler") or "").strip()
                elif field == "note":
                    new_row[field] = f"scheduler: {scheduler_note}"
                else:
                    # counts or any other extra columns
                    new_row[field] = ""

            summary_rows.append(new_row)
            index[(sra, srr)] = len(summary_rows) - 1
            continue

        # Existing row: append scheduler information to note.
        existing = (summary_rows[row_idx].get("note") or "").strip()
        if existing:
            merged = f"{existing} | scheduler: {scheduler_note}"
        else:
            merged = f"scheduler: {scheduler_note}"

        summary_rows[row_idx]["note"] = merged


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #
def run(summary_path: Path, trace_path: Path) -> int:
    """Core driver: load, annotate, write summary.tsv."""
    if not summary_path.is_file():
        LOG.error("Summary file not found: %s", summary_path)
        return 1

    if not trace_path.is_file():
        LOG.error("Trace file not found: %s", trace_path)
        return 1

    summary_rows, fieldnames = load_summary(summary_path)
    reasons = collect_reasons_from_trace(trace_path)

    if not reasons:
        LOG.info("No scheduler annotations to add; exiting")
        return 0

    metadata_root = summary_path.parent / "metadata"
    if not metadata_root.exists():
        LOG.warning(
            "Metadata directory %s does not exist; new rows will be "
            "created without platform/model/strategy/assembler fields",
            metadata_root,
        )

    annotate_summary(summary_rows, fieldnames, reasons, metadata_root)
    write_summary(summary_path, summary_rows, fieldnames)
    return 0


def main(argv: list[str] | None = None) -> int:
    """Entry point: parse args, configure logging, run pipeline."""
    args = parse_args(argv)
    setup_logging(args.verbose)

    summary_path = args.summary_tsv.expanduser().resolve()
    trace_path = args.trace_txt.expanduser().resolve()

    LOG.info("Summary file: %s", summary_path)
    LOG.info("Trace file:   %s", trace_path)

    return run(summary_path, trace_path)


if __name__ == "__main__":
    sys.exit(main())
