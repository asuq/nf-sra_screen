#!/usr/bin/env bash
# Run Rosella on an assembly + BAM
#
# Args:
#   --assembly     assembly fasta
#   --bam          BAM mapped against the assembly
#   --cpus         threads to use
#   --attempt      current attempt number (for Nextflow retries)
#   --max-retries  maximum attempts before treating failures as soft and writing rosella.note

set -euo pipefail

export MPLCONFIGDIR="$PWD/.mplconfig"
export NUMBA_CACHE_DIR="$PWD/.numba_cache"

mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"

assembly=""
bam=""
cpus=1
attempt=0
max_retries=1
has_bins=false 

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

tmp_dir="tmp_rosella"
final_dir="rosella"
log_file="${tmp_dir}/rosella_recover.log"
note_file="rosella.note"

rm -rf "$tmp_dir" "$final_dir" "$note_file"
mkdir -p "$tmp_dir" "$final_dir"
: > "$note_file"

fail() {
  local msg="$1"
  echo "$msg" >&2
  if [[ "$attempt" -lt "$max_retries" ]]; then
    exit 1
  fi
  rm -rf "$final_dir"
  mkdir -p "$final_dir"
  printf '%s\n' "$msg" > "$note_file"
  exit 0
}


if ! rosella recover \
        --assembly "$assembly" \
        --output-directory "$tmp_dir" \
        --bam-files "$bam" \
        --threads "$cpus" \
        &> "$log_file"; then

  # Special-case: known Rosella bug when there are 0 bins to refine
  if grep -q "attempt to divide by zero" "$log_file"; then
    echo "Rosella: no bins were generated" >&2
    mkdir -p "$final_dir"
    exit 0
  fi

  fail "Rosella: recover failed"
fi

mkdir -p "${tmp_dir}/unbinned"

shopt -s nullglob

# Move unbinned sequences aside
for f in "$tmp_dir"/*unbinned.fna; do
  [[ -e "$f" ]] || break
  mv "$f" "${tmp_dir}/unbinned/"
done

# Move bins if any
for f in "$tmp_dir"/rosella_*.fna; do
  [[ -e "$f" ]] || break
  has_bins=true
  mv "$f" "$final_dir"/
done

shopt -u nullglob

if [[ "$has_bins" = false ]]; then
  echo "Rosella: no rosella_*.fna bins produced" >&2
fi
