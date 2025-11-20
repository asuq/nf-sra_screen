#!/usr/bin/env bash
# Usage:
#   run_metaflye_nano.sh \
#     --reads READS \
#     --cpus N \
#     --attempt A \
#     --max-retries M
#
# Produces:
#   - assembly.fasta
#   - assembly.bam
#   - assembly.bam.csi
#   - assembly.gfa
#		- flye.log
#   - FAIL.note (only on fatal problems)

set -euo pipefail

read_files=()
cpus=1
attempt=0
max_retries=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --reads)
      shift
      while [[ $# -gt 0 ]]; do
        read_files+=( "$1" )
        shift
      done
      ;;
    --cpus) cpus="$2"; shift 2 ;;
    --attempt) attempt="$2"; shift 2 ;;
    --max-retries) max_retries="$2"; shift 2 ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1 ;;
  esac
done

if [[ ${#read_files[@]} -eq 0 ]]; then
  echo "run_metaflye_nano.sh: missing --reads" >&2
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

if [[ -z "${reads}" ]]; then
  echo "run_metaflye_nano.sh: missing --reads" >&2
  exit 1
fi

# Run metaFlye (ONT)
if ! flye --nano-raw "${read_files[@]}" \
          --threads "${cpus}" --scaffold --out-dir '.' --meta; then
  fail "metaFlye (ONT): assembly failed"
fi

# Map reads back with minimap2 + samtools
if ! ( minimap2 -ax map-ont -t "${cpus}" assembly.fasta "${read_files[@]}" \
      | samtools sort --output-fmt BAM -@ "${cpus}" -o assembly.bam \
      && samtools index -c -o assembly.bam.csi -@ "${cpus}" assembly.bam ); then
  fail "metaFlye (ONT): mapping/indexing failed"
fi

# Rename graph to standard name
mv -v assembly_graph.gfa assembly.gfa
