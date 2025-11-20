#!/usr/bin/env bash
# sandpiper_lookup.sh
# Usage: sandpiper_lookup.sh <sra_accession> <sandpiper_db_ch>
#  - <sra_accession>: SRR/ERR/etc. accession ID
#  - <sandpiper_db_ch>: directory containing sandpiper_sra.txt and sandpiper1.0.0.condensed.tsv

set -euo pipefail

if [ "$#" -ne 2 ]; then
    printf 'Usage: %s <sra_accession> <sandpiper_db_ch>\n' "$0" >&2
    exit 1
fi

ACCESSION=$1
DB_DIR=$2

# Remove a possible trailing slash from DB_DIR so we don't get // in paths
DB_DIR_STRIPPED=${DB_DIR%/}

SRA_LIST_FILE="${DB_DIR_STRIPPED}/sandpiper_sra.txt"
SANDPIPER_TSV="${DB_DIR_STRIPPED}/sandpiper1.0.0.condensed.tsv"

# Check that required files exist
if [ ! -f "$SRA_LIST_FILE" ]; then
    printf 'Error: SRA list file "%s" not found.\n' "$SRA_LIST_FILE" >&2
    exit 1
fi

if [ ! -f "$SANDPIPER_TSV" ]; then
    printf 'Error: sandpiper TSV file "%s" not found.\n' "$SANDPIPER_TSV" >&2
    exit 1
fi

# Use sandpiper_sra.txt to decide if this accession has a sandpiper result
# ^ACC$ anchors to the whole line; safe for SRA-style accessions
if ! grep "^${ACCESSION}\$" "$SRA_LIST_FILE" >/dev/null 2>&1; then
    printf 'no sandpiper result\n'
    exit 0
fi

# Accession is present in the SRA list: extract matching line(s) from the TSV
printf 'sample\tcoverage\ttaxonomy\n'
awk -F '\t' -v acc="$ACCESSION" '
    $1 == acc { print }
' "$SANDPIPER_TSV"
exit 0
