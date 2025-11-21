#!/usr/bin/env bash
# Usage:
#   run_download_srr.sh \
#     --srr SRR \
#     --platform PLATFORM \
#     --assembler ASM \
#     --cpus N \
#     --attempt A \
#     --max-retries M
#
# Produces:
#   - *.fastq.gz
#   - assembler.txt
#   - FAIL.note (only on fatal problems)

set -euo pipefail

srr=""
platform=""
assembler=""
cpus=1
attempt=0
max_retries=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --srr) srr="$2"; shift 2 ;;
    --platform) platform="$2"; shift 2 ;;
    --assembler) assembler="$2"; shift 2 ;;
    --cpus) cpus="$2"; shift 2 ;;
    --attempt) attempt="$2"; shift 2 ;;
    --max-retries) max_retries="$2"; shift 2 ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1 ;;
  esac
done

if [[ -z "${srr}" || -z "${platform}" ]]; then
  echo "run_download_srr.sh: missing required arguments (--srr, --platform)" >&2
  exit 1
fi

# Download sequence data
if ! iseq -Q -g -t "${cpus}" -p 8 -i "${srr}"; then
  if [[ "${attempt}" -lt "${max_retries}" ]]; then
    # Ask Nextflow to retry
    exit 1
  fi
  echo "Fastq: download raw data failed" > FAIL.note
  # Clean up partial SRA / FASTQ
  rm -f ./*.f*q* "${srr}" 2>/dev/null || true
  exit 0
fi

# PacBio assembler check (only if platform is PACBIO_SMRT and assembler is missing/unknown)
final_asm="${assembler}"
if [[ "${platform}" == "PACBIO_SMRT" && ( -z "${assembler}" || "${assembler}" == "unknown" ) ]]; then
  echo "Checking PacBio reads to determine assembler" >&2
  if zcat -f ./*.f*q* 2>/dev/null \
    | awk 'NR%4==1{ h=tolower($0); if (h ~ /\/ccs([[:space:]]|$)/) { found=1; exit } } END{ exit(!found) }'
  then
    final_asm="long_hifi"
  else
    final_asm="long_pacbio"
  fi
  printf '%s\n' "${final_asm}" > assembler.txt
else
  # Always provide an assembler.txt for downstream mapping
  : > assembler.txt
fi
