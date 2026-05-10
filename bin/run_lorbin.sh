#!/usr/bin/env bash
# Run LorBin binning on an assembly + BAM.
#
# Args:
#   --assembly    assembly fasta
#   --bam         BAM mapped against the assembly
#   --cpus        threads to use
#   --attempt     current attempt number (for Nextflow retries)
#   --max-retries maximum attempts before treating failures as soft
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
      echo "run_lorbin.sh: unknown argument: $1" >&2
      exit 1 ;;
  esac
done

if [[ -z "$assembly" || -z "$bam" ]]; then
  echo "run_lorbin.sh: missing --assembly or --bam" >&2
  exit 1
fi

tmp_dir="tmp_lorbin"
final_dir="lorbin"
note_file="lorbin.note"
map_file="lorbin.contig2bin.tsv"

rm -rf "$tmp_dir" "$final_dir" "$note_file" "$map_file"
: > "$note_file"
: > "$map_file"

record_soft_failure() {
  local msg="$1"
  rm -rf "$final_dir"
  mkdir -p "$final_dir"
  : > "$map_file"
  printf '%s\n' "$msg" > "$note_file"
}

handle_unexpected_exit() {
  local exit_code=$?

  if (( exit_code == 0 || attempt <= max_retries )); then
    return
  fi

  record_soft_failure "LorBin: unexpected failure (exit ${exit_code})"
  exit 0
}

trap handle_unexpected_exit EXIT

fail() {
  local msg="$1"
  echo "$msg" >&2
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
        print("LorBin GPU mode requires PyTorch, but torch is not installed", file=sys.stderr)
        return 1

    import torch

    if not torch.cuda.is_available():
        print("LorBin GPU mode requested, but CUDA is not visible to PyTorch", file=sys.stderr)
        return 1

    print(
        f"LorBin GPU preflight: CUDA visible to PyTorch "
        f"(cuda={torch.version.cuda}, devices={torch.cuda.device_count()})",
        file=sys.stderr,
    )
    return 0


raise SystemExit(main())
PY
}

if [[ "$require_cuda" == true ]]; then
  require_pytorch_cuda || fail "LorBin: GPU mode requested but CUDA preflight failed"
fi

mkdir -p "$tmp_dir"

if ! python - "$tmp_dir" "$assembly" "$bam" "$cpus" <<'PY'; then
"""Run LorBin with numeric thread arguments restored."""
import sys

from lorbin import lorbin


def parser_args_with_numeric_threads(original_parser_args):
    """Return parsed LorBin arguments after casting numeric CLI values."""
    args = original_parser_args()
    args.num_process = int(args.num_process)
    args.bin_length = int(args.bin_length)
    args.akeep = float(args.akeep)
    if not hasattr(args, "epoch"):
        args.epoch = 300
    return args


def patched_parser_args():
    """Parse LorBin CLI arguments with numeric fields fixed."""
    return parser_args_with_numeric_threads(original_parser_args)


output, fasta, bam, cpus = sys.argv[1:5]
original_parser_args = lorbin.parser_args
lorbin.parser_args = patched_parser_args
sys.argv = [
    "LorBin",
    "bin",
    "--output",
    output,
    "--fasta",
    fasta,
    "--bam",
    bam,
    "--num_process",
    str(int(cpus)),
]
lorbin.main()
PY
  fail "LorBin: bin command failed"
fi

bins_src="${tmp_dir}/output_bins"
if [[ ! -d "$bins_src" ]]; then
  fail "LorBin: output_bins directory not found"
fi

mkdir -p "$final_dir"
shopt -s nullglob
for f in "$bins_src"/*; do
  [[ -f "$f" ]] || continue
  cp "$f" "$final_dir"/
  awk -v bin="${f##*/}" '
    /^>/ {
      sub(/^>/, "")
      split($0, fields, /[[:space:]]+/)
      print fields[1] "\t" bin
    }
  ' "$f" >> "$map_file"
done
shopt -u nullglob
