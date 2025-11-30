#!/usr/bin/env bash
# Run Rosella on an assembly + BAM

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
      echo "run_rosella.sh: unknown argument: $1" >&2
      exit 1 ;;
  esac
done

if [[ -z "$assembly" || -z "$bam" ]]; then
  echo "run_rosella.sh: missing --assembly or --bam" >&2
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

tmp_dir="tmp_rosella"
final_dir="rosella"

rm -rf "$tmp_dir" "$final_dir"
mkdir -p "$tmp_dir"

if ! rosella recover \
        --assembly "$assembly" \
        --output-directory "$tmp_dir" \
        --bam-files "$bam" \
        --threads "$cpus"; then
  fail "Rosella: recover failed"
fi

mkdir -p "${tmp_dir}/unbinned" "$final_dir"

shopt -s nullglob

# Move unbinned sequences aside
for f in "$tmp_dir"/*unbinned.fna; do
  [[ -e "$f" ]] || break
  mv "$f" "${tmp_dir}/unbinned/"
done

# Move bins if any
has_bins=false
for f in "$tmp_dir"/rosella_*.fna; do
  [[ -e "$f" ]] || break
  has_bins=true
  mv "$f" "$final_dir"/
done

if [[ "$has_bins" = false ]]; then
  echo "Rosella: no rosella_*.fna bins produced (this is allowed)" >&2
fi

shopt -u nullglob
