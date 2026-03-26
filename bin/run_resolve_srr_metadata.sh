#!/usr/bin/env bash
# Usage:
#   run_resolve_srr_metadata.sh \
#     --sample SAMPLE \
#     --srr SRR \
#     --attempt A \
#     --max-retries M
#
# Produces:
#   - resolved.tsv
#   - FAIL.note (only on final fatal problems)

set -euo pipefail

sample=""
srr=""
attempt=0
max_retries=1

usage() {
  echo "Usage: $0 --sample SAMPLE --srr SRR --attempt A --max-retries M" >&2
}

fail() {
  local msg="$1"

  echo "$msg" >&2

  if [[ "${attempt}" -lt "${max_retries}" ]]; then
    exit 1
  fi

  printf '%s\n' "${msg}" > FAIL.note
  exit 0
}

parse_args() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --sample)
        sample="$2"
        shift 2
        ;;
      --srr)
        srr="$2"
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
        usage
        echo "Unknown argument: $1" >&2
        exit 1
        ;;
    esac
  done
}

run_shared_metadata() {
  if ! run_download_metadata.sh \
      --sra "${srr}" \
      --attempt "${attempt}" \
      --max-retries "${max_retries}"
  then
    exit 1
  fi

  if [[ -f FAIL.note ]]; then
    exit 0
  fi
}

write_resolved_metadata() {
  local filtered_csv="${srr}.filtered.csv"
  local skipped_csv="${srr}.skipped.csv"
  local line=""

  if [[ -s "${filtered_csv}" ]] && [[ "$(wc -l < "${filtered_csv}")" -ge 2 ]]; then
    line="$(awk -F',' -v run="${srr}" 'NR > 1 && $2 == run {print $3 "\t" $4 "\t" $6 "\t" $7; exit}' "${filtered_csv}")"
    if [[ -n "${line}" ]]; then
      printf '%s\n' "${line}" > resolved.tsv
    fi

    if [[ ! -s resolved.tsv ]]; then
      fail "Metadata: failed to resolve platform and assembler for sample '${sample}' run '${srr}'"
    fi

    return 0
  fi

  local skip_reason=""
  if [[ -s "${skipped_csv}" ]] && [[ "$(wc -l < "${skipped_csv}")" -ge 2 ]]; then
    skip_reason="$(awk -F',' -v run="${srr}" 'NR > 1 && $2 == run {print $7; exit}' "${skipped_csv}")"
  fi

  if [[ -n "${skip_reason}" ]]; then
    fail "Metadata: sample '${sample}' run '${srr}' was filtered out (${skip_reason})"
  fi

  fail "Metadata: no usable metadata remained for sample '${sample}' run '${srr}'"
}

main() {
  parse_args "$@"

  if [[ -z "${sample}" || -z "${srr}" ]]; then
    usage
    fail "Metadata: missing required arguments (--sample, --srr)"
  fi

  run_shared_metadata
  write_resolved_metadata
}

main "$@"
