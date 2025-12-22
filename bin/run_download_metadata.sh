#!/usr/bin/env bash
# Usage:
#   run_download_metadata.sh \
#       --sra ACCESSION \
#       --attempt A \
#       --max-retries M
#
# Produces:
#   - ${sra}.metadata.tsv/.csv (by iseq; used by filter_sra.sh)
#   - ${sra}.filtered.csv
#   - ${sra}.skipped.csv
#   - FAIL.note (when failed on final attempt)

set -euo pipefail

sra=""
attempt=0
max_retries=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sra)
      sra="$2"
      shift 2
      ;;
    --attempt)
      attempt="$2"
      shift 2
      ;;
    --max-retries)
      max_retries="$2"
      shift 2
      ;;
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
    # Let Nextflow retry
    exit 1
  fi
  echo "$msg" > FAIL.note
  exit 0
}

if [[ -z "${sra}" ]]; then
  fail "Metadata: missing --sra"
fi

# Retrieve metadata
if ! iseq -m -i "${sra}"; then
  fail "Metadata: iseq -m failed for ${sra}"
fi

# No metadata file
if [[ ! -f "${sra}.metadata.tsv" ]] && [[ ! -f "${sra}.metadata.csv" ]]; then
  fail "Metadata: no metadata file found for ${sra}"
fi

# Filter SRRs into .filtered.csv / .skipped.csv
if ! filter_sra.sh "${sra}"; then
  fail "Metadata: filter_sra.sh failed for ${sra}"
fi
