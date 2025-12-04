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

rm -rf "$tmp_dir" "$final_dir" "$note_file"
mkdir -p "$tmp_dir"
: > "$note_file"


fail() {
  local msg="$1"
  echo "$msg" >&2
  if [[ "$attempt" -lt "$max_retries" ]]; then
    exit 1
  fi
  # Final attempt: record soft fail but still emit an empty metabat dir & note
  rm -rf "$final_dir"
  mkdir -p "$final_dir"
  printf '%s\n' "$msg" > "$note_file"
  exit 0
}


if ! jgi_summarize_bam_contig_depths --outputDepth "depth.txt" "$bam"; then
  fail "Metabat: jgi_summarize_bam_contig_depths failed"
fi

if ! metabat2 --inFile "$assembly" \
              --abdFile "depth.txt" \
              --outFile "${tmp_dir}/metabatbin" \
              --numThreads "$cpus"; then
  fail "Metabat: metabat2 failed"
fi

mv "$tmp_dir" "$final_dir"
