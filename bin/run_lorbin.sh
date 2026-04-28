#!/usr/bin/env bash
# Run LorBin binning on an assembly + BAM.
#
# Args:
#   --assembly    assembly fasta
#   --bam         BAM mapped against the assembly
#   --cpus        threads to use
#   --attempt     current attempt number (for Nextflow retries)
#   --max-retries maximum attempts before treating failures as soft

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

fail() {
  local msg="$1"
  echo "$msg" >&2
  if [[ "$attempt" -lt "$max_retries" ]]; then
    exit 1
  fi
  rm -rf "$final_dir"
  mkdir -p "$final_dir"
  : > "$map_file"
  printf '%s\n' "$msg" > "$note_file"
  exit 0
}

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
