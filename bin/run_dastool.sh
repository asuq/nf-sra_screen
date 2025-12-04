#!/usr/bin/env bash
# Run DAS_Tool to integrate bins from MetaBAT, CONCOCT, SemiBin, and Rosella
#
# Args:
#   --assembly     assembly fasta (required)
#   --metabat-dir  MetaBAT bins directory (may be empty or non-existent)
#   --concoct-dir  CONCOCT bins directory (may be empty or non-existent)
#   --semibin-dir  SemiBin bins directory (unused but kept for symmetry)
#   --semibin-map  SemiBin contig_bins.tsv (contig-to-bin mapping)
#   --rosella-dir  Rosella bins directory (may be empty or non-existent)
#   --cpus         threads to use
#   --attempt      current attempt number (for Nextflow retries)
#   --max-retries  maximum attempts before treating failures as soft and writing dastool.note

set -euo pipefail

assembly=""
metabat_dir=""
concoct_dir=""
semibin_dir=""
semibin_map=""
rosella_dir=""
cpus=1
attempt=0
max_retries=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --assembly)    assembly="$2";    shift 2 ;;
    --metabat-dir) metabat_dir="$2"; shift 2 ;;
    --concoct-dir) concoct_dir="$2"; shift 2 ;;
    --semibin-dir) semibin_dir="$2"; shift 2 ;;
    --semibin-map) semibin_map="$2"; shift 2 ;;
    --rosella-dir) rosella_dir="$2"; shift 2 ;;
    --cpus)        cpus="$2";        shift 2 ;;
    --attempt)     attempt="$2";     shift 2 ;;
    --max-retries) max_retries="$2"; shift 2 ;;
    *)
      echo "run_dastool.sh: unknown argument: $1" >&2
      exit 1 ;;
  esac
done

if [[ -z "$assembly" ]]; then
  echo "run_dastool.sh: --assembly is required" >&2
  exit 1
fi

tmp_dir="tmp_dastool"
final_dir="dastool"
note_file="dastool.note"

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

shopt -s nullglob

# Build contig2bin TSVs *only* for tools that actually have outputs
metabat_tsv="${tmp_dir}/metabat_bins.tsv"
concoct_tsv="${tmp_dir}/concoct_bins.tsv"
rosella_tsv="${tmp_dir}/rosella_bins.tsv"
semibin_tsv="${tmp_dir}/semibin_bins.tsv"

: > "$metabat_tsv"
: > "$concoct_tsv"
: > "$rosella_tsv"
: > "$semibin_tsv"

# MetaBAT: .fa bins in $metabat_dir
if [[ -n "$metabat_dir" && -d "$metabat_dir" ]]; then
  for f in "$metabat_dir"/*.fa; do
    [[ -e "$f" ]] || break
    # contig_id<TAB>bin_filename
    grep '^>' "$f" \
      | cut -f1 \
      | sed 's/^>//' \
      | sed "s/\$/\t${f##*/} /" \
      >> "$metabat_tsv"
  done
fi

# CONCOCT: .fa bins in $concoct_dir
if [[ -n "$concoct_dir" && -d "$concoct_dir" ]]; then
  for f in "$concoct_dir"/*.fa; do
    [[ -e "$f" ]] || break
    grep '^>' "$f" \
      | sed "s/>.*/&\t${f##*/} /" \
      | sed 's/^>//' \
      >> "$concoct_tsv"
  done
fi

# Rosella: rosella_*.fna bins in $rosella_dir
if [[ -n "$rosella_dir" && -d "$rosella_dir" ]]; then
  for f in "$rosella_dir"/rosella_*.fna; do
    [[ -e "$f" ]] || break
    grep '^>' "$f" \
      | sed "s/>.*/&\t${f##*/} /" \
      | sed 's/^>//' \
      >> "$rosella_tsv"
  done
fi

# SemiBin: contig_bins.tsv (header + 2 columns)
if [[ -n "$semibin_map" && -f "$semibin_map" ]]; then
  # Drop header and reuse as the contig2bin mapping
  tail -n +2 "$semibin_map" >> "$semibin_tsv"
fi

# If no TSV has any content, we consider this a "no bins" situation.
bins_list=()
[[ -s "$metabat_tsv"  ]] && bins_list+=("$metabat_tsv")
[[ -s "$concoct_tsv"  ]] && bins_list+=("$concoct_tsv")
[[ -s "$rosella_tsv"  ]] && bins_list+=("$rosella_tsv")
[[ -s "$semibin_tsv"  ]] && bins_list+=("$semibin_tsv")

if (( ${#bins_list[@]} == 0 )); then
  msg="DASTool: no bins found for any tool"
  echo "$msg" >&2
  printf '%s\n' "$msg" > "$note_file"
  exit 0
fi

# Build comma‑separated list for DAS_Tool
bins_arg=$(IFS=,; echo "${bins_list[*]}")

dastool_log="${tmp_dir}/dastool.log"
if ! DAS_Tool \
      --bins "$bins_arg" \
      --contigs "$assembly" \
      --outputbasename "${tmp_dir}/dastool" \
      --write_bin_evals \
      --write_bins \
      --threads "$cpus" \
      >"$dastool_log" 2>&1; then

  # Special case: no high-quality bins – not a real error for our pipeline
  if grep -q "No bins with bin-score >0.5 found" "$dastool_log"; then
    msg="DASTool: no high-quality bins (bin-score > 0.5)"
    echo "$msg" >&2
    printf '%s\n' "$msg" > "$note_file"
    exit 0
  fi

  # Any other DASTool failure is treated as a genuine error
  fail "DASTool: refinement failed"
fi

if [[ ! -d "${tmp_dir}/dastool_DASTool_bins" ]]; then
  fail "DASTool: dastool_DASTool_bins directory not found"
fi

# Move resulting bins to final directory
for f in "${tmp_dir}/dastool_DASTool_bins"/*; do
  [[ -e "$f" ]] || break
  mv "$f" "$final_dir"/
done
shopt -u nullglob
