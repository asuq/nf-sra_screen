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

write_metadata_tsv() {
  local metadata_tsv="${sample}.metadata.tsv"
  local url="https://www.ebi.ac.uk/ena/portal/api/filereport"
  local fields="run_accession,instrument_platform,instrument_model,library_source,library_strategy"

  curl -fsSL \
    --retry 3 \
    --retry-delay 2 \
    --connect-timeout 20 \
    --max-time 60 \
    "${url}?accession=${srr}&result=read_run&fields=${fields}&format=tsv" \
    > "${metadata_tsv}"

  if [[ ! -s "${metadata_tsv}" ]]; then
    fail "Metadata: ENA returned no metadata for ${srr}"
  fi

  local line_count
  line_count="$(wc -l < "${metadata_tsv}")"
  if [[ "${line_count}" -lt 2 ]]; then
    fail "Metadata: ENA returned no run rows for ${srr}"
  fi

  local resolved_srr
  resolved_srr="$(awk -F'\t' 'NR==2 {print $1}' "${metadata_tsv}")"
  if [[ -z "${resolved_srr}" ]]; then
    fail "Metadata: missing run_accession for ${srr}"
  fi

  if [[ "${resolved_srr}" != "${srr}" ]]; then
    fail "Metadata: resolved run_accession '${resolved_srr}' does not match requested SRR '${srr}'"
  fi
}

write_resolved_metadata() {
  local filtered_csv="${sample}.filtered.csv"
  local skipped_csv="${sample}.skipped.csv"

  if ! filter_sra.sh "${sample}"; then
    fail "Metadata: filter_sra.sh failed for ${srr}"
  fi

  if [[ -s "${filtered_csv}" ]] && [[ "$(wc -l < "${filtered_csv}")" -ge 2 ]]; then
    awk -F',' 'NR==2 {print $3 "\t" $4 "\t" $6 "\t" $7}' "${filtered_csv}" > resolved.tsv

    if [[ ! -s resolved.tsv ]]; then
      fail "Metadata: failed to resolve platform and assembler for ${srr}"
    fi

    return 0
  fi

  local skip_reason=""
  if [[ -s "${skipped_csv}" ]] && [[ "$(wc -l < "${skipped_csv}")" -ge 2 ]]; then
    skip_reason="$(awk -F',' 'NR==2 {print $7}' "${skipped_csv}")"
  fi

  if [[ -n "${skip_reason}" ]]; then
    fail "Metadata: SRR ${srr} was filtered out (${skip_reason})"
  fi

  fail "Metadata: no usable metadata remained for ${srr}"
}

main() {
  parse_args "$@"

  if [[ -z "${sample}" || -z "${srr}" ]]; then
    usage
    fail "Metadata: missing required arguments (--sample, --srr)"
  fi

  write_metadata_tsv
  write_resolved_metadata
}

main "$@"
