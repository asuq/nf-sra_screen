#!/usr/bin/env bash
# Run MetaBAT2 binning on an assembly + BAM
#
# Args:
#   --assembly   assembly fasta
#   --bam        BAM mapped against the assembly
#   --cpus       threads to use
#   --attempt    current attempt number (for Nextflow retries)
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

fail() {
  local msg="$1"
  echo "$msg" >&2
  if [[ "$attempt" -lt "$max_retries" ]]; then
    exit 1
  fi
  echo "$msg" > FAIL.note
  exit 0
}

tmp_dir="tmp_metabat"
out_dir="metabat_tmp"
final_dir="metabat"

rm -rf "$tmp_dir" "$out_dir" "$final_dir"
mkdir -p "$tmp_dir" "$out_dir"

if ! jgi_summarize_bam_contig_depths --outputDepth "${tmp_dir}/depth.txt" "$bam"; then
  fail "Metabat: jgi_summarize_bam_contig_depths failed"
fi

if ! metabat2 --inFile "$assembly" \
              --abdFile "${tmp_dir}/depth.txt" \
              --outFile "${out_dir}/metabatbin" \
              --numThreads "$cpus"; then
  fail "Metabat: metabat2 failed"
fi

mv "$out_dir" "$final_dir"
