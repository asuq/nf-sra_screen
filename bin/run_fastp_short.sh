#!/usr/bin/env bash
# Run fastp for assembler-scoped short-read preprocessing.

set -euo pipefail

srr=""
assembler=""
threads=1
attempt=0
max_retries=1
read_files=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --reads)
      shift
      while [[ $# -gt 0 ]]; do
        case "$1" in
          --srr|--assembler|--cpus|--attempt|--max-retries)
            break
            ;;
          *)
            read_files+=("$1")
            shift
            ;;
        esac
      done
      ;;
    --srr)         srr="$2"; shift 2 ;;
    --assembler)   assembler="$2"; shift 2 ;;
    --cpus)        threads="$2"; shift 2 ;;
    --attempt)     attempt="$2"; shift 2 ;;
    --max-retries) max_retries="$2"; shift 2 ;;
    *)
      echo "run_fastp_short.sh: unknown argument: $1" >&2
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
  printf '%s\n' "$msg" > FAIL.note
  exit 0
}

if [[ -z "$srr" || -z "$assembler" ]]; then
  echo "run_fastp_short.sh: --srr and --assembler are required" >&2
  exit 1
fi

if (( ${#read_files[@]} < 2 )); then
  fail "fastp: paired-end reads not found"
fi

length_required=50
if [[ "$assembler" == "metaspades" ]]; then
  length_required=100
fi

if ! fastp \
    --in1 "${read_files[0]}" \
    --in2 "${read_files[1]}" \
    --out1 "${srr}_fastp_R1.fastq.gz" \
    --out2 "${srr}_fastp_R2.fastq.gz" \
    --thread "$threads" \
    --length_required "$length_required" \
    --detect_adapter_for_pe \
    --qualified_quality_phred 30 \
    --html fastp.html \
    --json fastp.json; then
  fail "fastp: preprocessing failed"
fi
