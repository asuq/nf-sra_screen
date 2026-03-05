#!/usr/bin/env bash
# Usage:
#   run_metaflye_hifi.sh \
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
#   - flye.log
#   - FAIL.note (only on fatal problems)

set -euo pipefail

cpus=1
attempt=0
max_retries=1
read_files=()

cleanup() {
  local exit_status=$?

  # Never let cleanup failures replace the real script status.
  set +e

  case "${work_dir:-}" in
    ""|"/")
      printf 'Refusing cleanup in unsafe directory: %q\n' "${work_dir:-<unset>}" >&2
      return "$exit_status"
      ;;
  esac

  find "$work_dir" \
    -mindepth 1 \
    -maxdepth 1 \
    -type d \
    -exec rm -rf -- {} +

  return "$exit_status"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --reads)
      shift
      while [[ $# -gt 0 ]]; do
        # Stop when we hit the next flag
        case "$1" in
          --sandpiper-decision|--valid-taxa|--singlem-db|--cpus|--attempt|--max-retries)
            break
            ;;
          *)
            read_files+=("$1")
            shift
            ;;
        esac
      done
      ;;
    --cpus)
      cpus="$2"
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

if [[ ${#read_files[@]} -eq 0 ]]; then
  echo "run_metaflye_hifi.sh: missing --reads" >&2
  exit 1
fi

readonly work_dir=$PWD
trap cleanup EXIT

fail() {
  local msg="$1"
  echo "$msg" >&2
  if [[ "$attempt" -lt "$max_retries" ]]; then
    exit 1
  fi
  echo "$msg" > FAIL.note
  exit 0
}

# Run metaFlye (HiFi)
if ! flye --pacbio-hifi "${read_files[@]}" \
          --threads "${cpus}" \
          --scaffold \
          --out-dir '.' \
          --meta; then
  fail "metaFlye (HiFi): assembly failed"
fi

# Map reads back with minimap2 + samtools
if ! (
  minimap2 -ax map-hifi -I 20G -t "${cpus}" assembly.fasta "${read_files[@]}" \
    | samtools sort --output-fmt BAM -@ "${cpus}" -o assembly.bam \
  && samtools index -c -o assembly.bam.csi -@ "${cpus}" assembly.bam
); then
  fail "metaFlye (HiFi): mapping/indexing failed"
fi

# Rename graph to standard name
if ! mv -v assembly_graph.gfa assembly.gfa; then
  fail "metaFlye (HiFi): graph rename failed"
fi