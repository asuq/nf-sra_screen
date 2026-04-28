#!/usr/bin/env bash
# Run SemiBin2 single_easy_bin on an assembly + BAM
# Uses DIAMONDDB from --diamond-db (same DB as DIAMOND process)
#
# Args:
#   --assembly     assembly fasta
#   --bam          BAM mapped against the assembly
#   --diamond-db   DIAMOND reference database
#   --cpus         threads to use
#   --attempt      current attempt number (for Nextflow retries)
#   --max-retries  maximum attempts before treating failures as soft and writing semibin.note
#   --require-cuda fail before binning if PyTorch cannot see CUDA

set -euo pipefail

assembly=""
bam=""
diamond_db=""
cpus=1
attempt=0
max_retries=1
require_cuda=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --assembly) assembly="$2"; shift 2 ;;
    --bam) bam="$2"; shift 2 ;;
    --diamond-db) diamond_db="$2"; shift 2 ;;
    --cpus) cpus="$2"; shift 2 ;;
    --attempt) attempt="$2"; shift 2 ;;
    --max-retries) max_retries="$2"; shift 2 ;;
    --require-cuda) require_cuda=true; shift ;;
    *)
      echo "run_semibin.sh: unknown argument: $1" >&2
      exit 1 ;;
  esac
done

if [[ -z "$assembly" || -z "$bam" || -z "$diamond_db" ]]; then
  echo "run_semibin.sh: missing --assembly, --bam, or --diamond-db" >&2
  exit 1
fi

export DIAMONDDB="$diamond_db"

tmp_dir="tmp_semibin"
final_dir="semibin"
note_file="semibin.note"
map_file="semibin.contig2bin.tsv"

rm -rf "$tmp_dir" "$final_dir" "$note_file" "$map_file"
mkdir -p "$tmp_dir" "$final_dir"
: > "$note_file"
: > "$map_file"

fail() {
  local msg="$1"
  echo "$msg" >&2
  if [[ "$attempt" -lt "$max_retries" ]]; then
    exit 1
  fi
  rm -rf "$final_dir"
  mkdir -p "$final_dir"
  touch "$final_dir/contig_bins.tsv"
  : > "$map_file"
  printf '%s\n' "$msg" > "$note_file"
  exit 0
}

require_pytorch_cuda() {
  python - <<'PY'
"""Fail clearly unless PyTorch can access a CUDA device."""
import importlib.util
import sys


def main() -> int:
    """Check PyTorch CUDA visibility."""
    if importlib.util.find_spec("torch") is None:
        print("SemiBin2 GPU mode requires PyTorch, but torch is not installed", file=sys.stderr)
        return 1

    import torch

    if not torch.cuda.is_available():
        print("SemiBin2 GPU mode requested, but CUDA is not visible to PyTorch", file=sys.stderr)
        return 1

    print(
        f"SemiBin2 GPU preflight: CUDA visible to PyTorch "
        f"(cuda={torch.version.cuda}, devices={torch.cuda.device_count()})",
        file=sys.stderr,
    )
    return 0


raise SystemExit(main())
PY
}

if [[ "$require_cuda" == true ]]; then
  require_pytorch_cuda || fail "SemiBin2: GPU mode requested but CUDA preflight failed"
fi

if ! SemiBin2 single_easy_bin \
      --input-fasta "$assembly" \
      --input-bam "$bam" \
      --environment global \
      --output "$tmp_dir" \
      --threads "$cpus"; then
  fail "SemiBin2: single_easy_bin failed"
fi

if [[ ! -f "${tmp_dir}/contig_bins.tsv" ]]; then
  fail "SemiBin2: contig_bins.tsv not found"
fi

cp "${tmp_dir}/contig_bins.tsv" "${final_dir}/contig_bins.tsv"
tail -n +2 "${tmp_dir}/contig_bins.tsv" > "$map_file"

shopt -s nullglob

# Collect any bins produced (gzipped or not)
bins=( "${tmp_dir}/output_bins"/*.gz "${tmp_dir}/output_bins"/* )

if [[ "${#bins[@]}" -ne 0 ]]; then
  # Copy bins to final directory
  for f in "${bins[@]}"; do
    cp "$f" "$final_dir"/
  done

  # Parallel gunzip, if needed
  gz_files=( "$final_dir"/*.gz )
  if [[ "${#gz_files[@]}" -gt 0 ]]; then
    printf '%s\0' "${gz_files[@]}" \
      | xargs -0 -n 1 -P "${cpus}" gunzip -v \
      || fail "SemiBin2: gunzip failed"
  fi
fi

shopt -u nullglob
