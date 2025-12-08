#!/usr/bin/env bash
# Run COMEBin binning on an assembly + BAM
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
note_file="comebin.note"

rm -rf "$tmp_out" "$final_dir" "$note_file" "$bam_dir"
mkdir -p "$tmp_out" "$bam_dir"
: > "$note_file"

fail() {
  local msg="$1"
  echo "$msg" >&2
  if [[ "$attempt" -lt "$max_retries" ]]; then
    exit 1
  fi
  # Final attempt: record soft fail but still emit an empty comebin dir & note
  rm -rf "$final_dir"
  mkdir -p "$final_dir"
  printf '%s\n' "$msg" > "$note_file"
  exit 0
}

# Symlink the BAM into a directory, as COMEBin expects a bam directory (-p)
ln -s "$bam" "${bam_dir}/$(basename "$bam")"

# Run COMEBin (CPU mode)
if ! run_comebin.sh \
      -a "$assembly" \
      -p "$bam_dir" \
      -o "$tmp_out" \
      -n 6 \
      -t "$cpus"
then
  fail "COMEBin: run_comebin.sh failed"
fi

bins_src="${tmp_out}/comebin_res/comebin_res_bins"

if [[ ! -d "$bins_src" ]]; then
  fail "COMEBin: expected bins directory '${bins_src}' not found"
fi

mv "$bins_src" "$final_dir"
