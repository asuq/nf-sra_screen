#!/usr/bin/env bash
# Run COMEBin binning on an assembly + BAM
#
# Args:
#   --assembly    assembly fasta
#   --bam         BAM mapped against the assembly
#   --cpus        threads to use
#   --attempt     current attempt number (for Nextflow retries)
#   --max-retries maximum attempts before writing FAIL.note and exiting 0
#   --require-cuda fail before binning if PyTorch cannot see CUDA

set -euo pipefail

assembly=""
bam=""
cpus=1
attempt=0
max_retries=1
require_cuda=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --assembly) assembly="$2"; shift 2 ;;
    --bam) bam="$2"; shift 2 ;;
    --cpus) cpus="$2"; shift 2 ;;
    --attempt) attempt="$2"; shift 2 ;;
    --max-retries) max_retries="$2"; shift 2 ;;
    --require-cuda) require_cuda=true; shift ;;
    *)
      echo "run_comebin.sh: unknown argument: $1" >&2
      exit 1 ;;
  esac
done

if [[ -z "$assembly" || -z "$bam" ]]; then
  echo "run_comebin.sh: missing --assembly or --bam" >&2
  exit 1
fi

tmp_out="tmp_comebin"
bam_dir="comebin_bam"
final_dir="comebin"
note_file="FAIL.note"
map_file="comebin.contig2bin.tsv"
min_contig_len=1000
max_batch_size=1024
filtered_assembly="${tmp_out}/assembly.comebin.min1001.fasta"
filter_stats="${tmp_out}/assembly.comebin.filter.tsv"

rm -rf "$tmp_out" "$final_dir" "$note_file" "$map_file" "$bam_dir"
mkdir -p "$tmp_out" "$bam_dir"
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

  record_soft_failure "COMEBin: unexpected failure (exit ${exit_code})"
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

require_pytorch_cuda() {
  python - <<'PY'
"""Fail clearly unless PyTorch can access a CUDA device."""
import importlib.util
import sys


def main() -> int:
    """Check PyTorch CUDA visibility."""
    if importlib.util.find_spec("torch") is None:
        print("COMEBin GPU mode requires PyTorch, but torch is not installed", file=sys.stderr)
        return 1

    import torch

    if not torch.cuda.is_available():
        print("COMEBin GPU mode requested, but CUDA is not visible to PyTorch", file=sys.stderr)
        return 1

    print(
        f"COMEBin GPU preflight: CUDA visible to PyTorch "
        f"(cuda={torch.version.cuda}, devices={torch.cuda.device_count()})",
        file=sys.stderr,
    )
    return 0


raise SystemExit(main())
PY
}

if [[ "$require_cuda" == true ]]; then
  if require_pytorch_cuda; then
    :
  else
    fail "COMEBin: GPU mode requested but CUDA preflight failed" "$?"
  fi
fi

if filter_comebin_fasta.py \
      --input "$assembly" \
      --output "$filtered_assembly" \
      --stats "$filter_stats" \
      --min-length "$min_contig_len"
then
  :
else
  fail "COMEBin: failed to filter assembly for contigs longer than ${min_contig_len} bp" "$?"
fi

total_contigs="$(awk '$1 == "total" { print $2 }' "$filter_stats")"
eligible_contigs="$(awk '$1 == "kept" { print $2 }' "$filter_stats")"

if [[ -z "$total_contigs" || -z "$eligible_contigs" ]]; then
  fail "COMEBin: failed to read filtered assembly counts"
fi

if (( eligible_contigs < 2 )); then
  record_soft_failure "COMEBin: skipped because only ${eligible_contigs} contig(s) are longer than ${min_contig_len} bp (total contigs: ${total_contigs}); at least 2 eligible contigs are required for training"
  exit 0
fi

batch_size="$eligible_contigs"
if (( batch_size > max_batch_size )); then
  batch_size="$max_batch_size"
fi

echo "COMEBin: using ${eligible_contigs}/${total_contigs} contigs longer than ${min_contig_len} bp with batch size ${batch_size}" >&2

# Symlink the BAM into a directory, as COMEBin expects a bam directory (-p)
ln -s "$(realpath "$bam")" "${bam_dir}/$(basename "$bam")"

# Run COMEBin (CPU mode)
if run_comebin.sh \
      -a "$filtered_assembly" \
      -p "$bam_dir" \
      -o "$tmp_out" \
      -n 6 \
      -t "$cpus" \
      -b "$batch_size"
then
  :
else
  fail "COMEBin: run_comebin.sh failed" "$?"
fi

bins_src="${tmp_out}/comebin_res/comebin_res_bins"

if [[ ! -d "$bins_src" ]]; then
  fail "COMEBin: expected bins directory '${bins_src}' not found"
fi

mv "$bins_src" "$final_dir"

shopt -s nullglob
for f in "$final_dir"/*; do
  [[ -f "$f" ]] || continue
  awk -v bin="${f##*/}" '
    /^>/ {
      sub(/^>/, "")
      split($0, fields, /[[:space:]]+/)
      print fields[1] "\t" bin
    }
  ' "$f" >> "$map_file"
done
shopt -u nullglob
