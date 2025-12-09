#!/usr/bin/env bash
# Usage:
#   run_extract_taxa.sh \
#       --blobtable blobtools.csv \
#       --fasta assembly.fasta \
#       --taxa validated_taxa.csv \
#       --taxdump /path/to/taxdump \
#       --attempt A \
#       --max-retries M
#
# Notes:
#   - --taxa is the validated taxa file from validate_taxa.py, with columns:
#       rank,taxa,ncbi_phylum,gtdb_phylum
#   - For extraction, only NCBI-style taxa (not GTDB-style) are usable, since
#     BlobTools/taxdump are NCBI-based.
#   - We derive an intermediate CSV 'ncbi_taxa_for_extraction.csv' with header
#     'rank,ncbi_taxa' and pass that to extract_records.py.
#   - If there are no NCBI taxa (i.e. all input taxa are GTDB style), we skip
#     extraction, emit a FAIL.note with a warning, and exit 0.

set -euo pipefail

blobtable=""
fasta=""
taxa=""
taxdump=""
attempt=0
max_retries=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --blobtable)   blobtable="$2"; shift 2 ;;
    --fasta)       fasta="$2"; shift 2 ;;
    --taxa)        taxa="$2"; shift 2 ;;
    --taxdump)     taxdump="$2"; shift 2 ;;
    --attempt)     attempt="$2"; shift 2 ;;
    --max-retries) max_retries="$2"; shift 2 ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

fail() {
  local msg="$1"
  echo "$msg" >&2
  if [[ "$attempt" -lt "$max_retries" ]]; then
    exit 1
  fi
  echo "$msg" > FAIL.note
  exit 0
}

skip_extraction() {
  local msg="$1"
  echo "$msg" >&2
  echo "$msg" > FAIL.note
  exit 0
}

# Basic argument check
if [[ -z "$blobtable" || -z "$fasta" || -z "$taxa" || -z "$taxdump" ]]; then
  fail "Extract_records: missing required arguments (--blobtable, --fasta, --taxa, --taxdump)"
fi

if [[ ! -f "$taxa" ]]; then
  fail "Extract_records: taxa file not found: ${taxa}"
fi

# Derive NCBI taxa list from validated_taxa.csv:
#   - header: rank,taxa,ncbi_phylum,gtdb_phylum
#   - keep only rows where 'taxa' is NOT GTDB-style (no d__/p__/c__/o__/f__/g__/s__).
#   - output header: rank,ncbi_taxa
ncbi_taxa_csv="ncbi_taxa_for_extraction.csv"
rm -f "$ncbi_taxa_csv"

if ! awk -F',' '
  BEGIN {
    col_rank = 0
    col_taxa = 0
  }

  NR == 1 {
    # Parse header: detect rank and taxa columns (case-insensitive).
    for (i = 1; i <= NF; i++) {
      h = $i
      gsub(/^[ \t"]+|[ \t"]+$/, "", h)
      lh = tolower(h)
      if (lh == "rank") {
        col_rank = i
      } else if (lh == "taxa" || lh == "ncbi_taxa") {
        col_taxa = i
      }
    }
    if (!col_rank || !col_taxa) {
      print "ERROR: taxa file must have columns 'rank' and 'taxa'/'ncbi_taxa'" > "/dev/stderr"
      exit 1
    }

    print "rank,ncbi_taxa"
    next
  }

  NR > 1 {
    rank = $col_rank
    taxa = $col_taxa

    # Strip outer quotes and whitespace.
    gsub(/^"|"$/, "", rank)
    gsub(/^"|"$/, "", taxa)
    sub(/^[ \t]+/, "", rank); sub(/[ \t]+$/, "", rank)
    sub(/^[ \t]+/, "", taxa); sub(/[ \t]+$/, "", taxa)

    if (rank == "" || taxa == "") {
      next
    }

    # GTDB style detection (case-insensitive): d__/p__/c__/o__/f__/g__/s__
    lt = tolower(taxa)
    if (lt ~ /^(d__|p__|c__|o__|f__|g__|s__)/) {
      # GTDB-style taxa cannot be used for NCBI-based extraction; skip.
      next
    }

    key = rank "," taxa
    if (!(key in seen)) {
      seen[key] = 1
      print key
    }
  }
' "$taxa" > "$ncbi_taxa_csv"; then
  fail "Extract_records: failed to derive NCBI taxa list from ${taxa}"
fi

# If the derived file has only the header, there are no NCBI taxa.
line_count="$(wc -l < "$ncbi_taxa_csv" | tr -d '[:space:]')"
if [[ -z "$line_count" || "$line_count" -le 1 ]]; then
  skip_extraction "Extract_records: skipping extraction because taxa list contains only GTDB-style taxa (no NCBI taxa for BlobTools/taxdump)"
fi

# Normal extraction using the derived NCBI taxa list
if ! extract_records.py --blobtable "$blobtable" \
      --fasta "$fasta" --taxa "$ncbi_taxa_csv" --taxdump "$taxdump"; then
  fail "Extract_records: run failed"
fi
