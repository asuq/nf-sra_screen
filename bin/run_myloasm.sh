#!/usr/bin/env bash
# Usage:
#   run_myloasm.sh \
#     --read-type nanopore|hifi \
#     --reads READS \
#     --cpus N \
#     --attempt A \
#     --max-retries M
#
# Produces:
#   - assembly.fasta
#   - assembly.gfa
#   - myloasm_(DATETIME).log
#   - FAIL.note (only on fatal problems)

set -euo pipefail

cpus=1
attempt=0
max_retries=1
read_type=""
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
          --sandpiper-decision|--valid-taxa|--singlem-db|--read-type|--cpus|--attempt|--max-retries)
            break
            ;;
          *)
            read_files+=("$1")
            shift
            ;;
        esac
      done
      ;;
    --read-type)    read_type="$2"; shift 2 ;;
    --cpus)        cpus="$2"; shift 2 ;;
    --attempt)     attempt="$2"; shift 2 ;;
    --max-retries) max_retries="$2"; shift 2 ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

if [[ ${#read_files[@]} -eq 0 ]]; then
  echo "run_myloasm.sh: missing --reads" >&2
  exit 1
fi

mode_args=()
case "$read_type" in
  nanopore)
    # Nanopore R10 is Myloasm's default mode.
    ;;
  hifi)
    mode_args+=(--hifi)
    ;;
  "")
    echo "run_myloasm.sh: missing --read-type" >&2
    exit 1
    ;;
  *)
    echo "run_myloasm.sh: unsupported --read-type '$read_type'; expected nanopore or hifi" >&2
    exit 1
    ;;
esac

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

# Run myloasm
if ! myloasm "${read_files[@]}" -o . -t "${cpus}" "${mode_args[@]}"; then
  fail "myloasm: assembly failed"
fi

# Rename the polished, filtered assembly to the pipeline output name.
if ! mv -v -- assembly_primary.fa assembly.fasta; then
  fail "myloasm: failed to rename assembly_primary.fa to assembly.fasta"
fi

# Reconcile the pre-polish graph with the final assembly. By default,
# annotate-gfa drops per-read alignment records because their coordinates no
# longer match the polished sequences.
if ! mylotools annotate-gfa \
  --gfa final_contig_graph.gfa \
  --fasta assembly.fasta \
  --output assembly.gfa \
  || [[ ! -s assembly.gfa ]]; then
  rm -f -- assembly.gfa
  fail "myloasm: GFA annotation failed"
fi
