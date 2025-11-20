#!/usr/bin/env bash
# Usage:
#   run_extract_taxa.sh \
#       --blobtable blobtools.csv \
#       --fasta assembly.fasta \
#       --taxa taxa.csv \
#       --taxdump /path/to/taxdump \
#       --attempt A \
#       --max-retries M

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

if ! extract_records.py --blobtable "$blobtable" \
      --fasta "$fasta" --taxa "$taxa" --taxdump "$taxdump"; then
  fail "Extract_records: run failed"
fi
