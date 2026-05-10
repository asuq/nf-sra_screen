#!/usr/bin/env bash
# Run MetaBAT2 binning on an assembly + BAM
#
# Args:
#   --assembly    assembly fasta
#   --bam         BAM mapped against the assembly
#   --cpus        threads to use
#   --attempt     current attempt number (for Nextflow retries)
#   --max-retries maximum attempts before writing FAIL.note and exiting 0

set -euo pipefail

assembly=""
bam=""
cpus=1
attempt=0
max_retries=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --assembly) assembly="$2"; shift 2 ;;
    --bam) bam="$2"; shift 2 ;;
    --cpus) cpus="$2"; shift 2 ;;
    --attempt) attempt="$2"; shift 2 ;;
    --max-retries) max_retries="$2"; shift 2 ;;
    *)
      echo "run_metabat.sh: unknown argument: $1" >&2
      exit 1 ;;
  esac
done

if [[ -z "$assembly" || -z "$bam" ]]; then
  echo "run_metabat.sh: missing --assembly or --bam" >&2
  exit 1
fi

tmp_dir="tmp_metabat"
final_dir="metabat"
note_file="metabat.note"
map_file="metabat.contig2bin.tsv"

rm -rf "$tmp_dir" "$final_dir" "$note_file" "$map_file"
mkdir -p "$tmp_dir"
: > "$note_file"
: > "$map_file"


record_soft_failure() {
  local msg="$1"
  rm -rf "$final_dir"
  mkdir -p "$final_dir"
  : > "$map_file"
  printf '%s\n' "$msg" > "$note_file"
}

is_scheduler_failure_exit() {
  local exit_code="$1"
  case "$exit_code" in
    137|139|140|143) return 0 ;;
    *) return 1 ;;
  esac
}

handle_unexpected_exit() {
  local exit_code=$?

  if (( exit_code == 0 || attempt <= max_retries )) || is_scheduler_failure_exit "$exit_code"; then
    return
  fi

  record_soft_failure "Metabat: unexpected failure (exit ${exit_code})"
  exit 0
}

trap handle_unexpected_exit EXIT

fail() {
  local msg="$1"
  local exit_code="${2:-1}"
  echo "$msg" >&2
  if is_scheduler_failure_exit "$exit_code"; then
    exit "$exit_code"
  fi
  if (( attempt <= max_retries )); then
    exit 1
  fi
  record_soft_failure "$msg"
  exit 0
}


if jgi_summarize_bam_contig_depths --outputDepth "depth.txt" "$bam"; then
  :
else
  fail "Metabat: jgi_summarize_bam_contig_depths failed" "$?"
fi

if metabat2 --inFile "$assembly" \
              --abdFile "depth.txt" \
              --outFile "${tmp_dir}/metabatbin" \
              --numThreads "$cpus"; then
  :
else
  fail "Metabat: metabat2 failed" "$?"
fi

mv "$tmp_dir" "$final_dir"

shopt -s nullglob
for f in "$final_dir"/*.fa; do
  awk -v bin="${f##*/}" '
    /^>/ {
      sub(/^>/, "")
      split($0, fields, /[[:space:]]+/)
      print fields[1] "\t" bin
    }
  ' "$f" >> "$map_file"
done
shopt -u nullglob
