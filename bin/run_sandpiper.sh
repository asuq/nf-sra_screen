#!/usr/bin/env bash
# Usage:
#   run_sandpiper.sh \
#       --srr SRR \
#       --valid-taxa validated_taxa.csv \
#       --db /path/to/sandpiper_db \
#       --attempt A \
#       --max-retries M
#
# Produces:
#   - gtdb_taxa_to_check.txt
#   - sandpiper_report.txt
#   - sandpiper_output.tsv (optional)
#   - sandpiper_decision.txt
#   - FAIL.note (only on fatal problems / negative)

set -euo pipefail

attempt=0
max_retries=1
srr=""
valid_taxa=""
db_dir=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --srr) srr="$2"; shift 2 ;;
    --valid-taxa) valid_taxa="$2"; shift 2 ;;
    --db) db_dir="$2"; shift 2 ;;
    --attempt) attempt="$2"; shift 2 ;;
    --max-retries) max_retries="$2"; shift 2 ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

hard_fail() {
  local msg="$1"
  echo "$msg" >&2
  if [[ "$attempt" -lt "$max_retries" ]]; then
    exit 1
  fi
  # On final attempt, fall back to SingleM without writing FAIL.note
  # (to avoid duplicate notes per sample in the global summary)
  echo "$msg"
  echo "RUN_SINGLEM" > sandpiper_decision.txt
  exit 0
}

if [[ -z "$srr" || -z "$valid_taxa" || -z "$db_dir" ]]; then
  hard_fail "Sandpiper: missing required arguments (--srr, --valid-taxa, --db)"
fi


# Extract gtdb_taxa
make_gtdb_taxa_list_from_validated_taxa.sh "$valid_taxa" gtdb_taxa_to_check.txt \
  || hard_fail "Sandpiper: failed to parse validated_taxa"

# Lookup SRR in Sandpiper DB
sandpiper_lookup.sh "$srr" "$db_dir" > sandpiper_report.txt \
  || hard_fail "Sandpiper: lookup failed"


# If Sandpiper has no precomputed profile for this sample, fall back to SingleM
if grep -q '^no sandpiper result$' sandpiper_report.txt; then
  echo "Sandpiper: no precomputed result for ${srr}; falling back to SingleM" >&2
  echo "RUN_SINGLEM" > sandpiper_decision.txt
  exit 0
fi

# Check taxa
rc=0
if check_singlem_taxa.py \
	-i sandpiper_report.txt \
	-g gtdb_taxa_to_check.txt \
	-o sandpiper_output.tsv; then
  rc=0
else
  rc=$?
fi

case "$rc" in
  0)
    echo "PASS" > sandpiper_decision.txt
    exit 0
    ;;

  1)
    hard_fail "Sandpiper: phylum check internal error"
    ;;

  2)
    echo "NEGATIVE" > sandpiper_decision.txt
    echo "Sandpiper: No target gtdb_taxa detected" > FAIL.note
    exit 0
    ;;

  *)
    hard_fail "SingleM phylum check returned unexpected status ${rc}"
    ;;
esac
