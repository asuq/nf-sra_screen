#!/usr/bin/env bash
# Run DAS_Tool to integrate bins from MetaBAT, COMEBin, SemiBin, and Rosella
#
# Args:
#   --assembly     assembly fasta (required)
#   --bin-maps     comma-separated contig-to-bin TSVs from selected binners
#   --metabat-dir  legacy MetaBAT bins directory (may be empty or non-existent)
#   --comebin-dir  legacy COMEBin bins directory (may be empty or non-existent)
#   --semibin-dir  legacy SemiBin bins directory (unused but kept for symmetry)
#   --semibin-map  legacy SemiBin contig_bins.tsv (contig-to-bin mapping)
#   --rosella-dir  legacy Rosella bins directory (may be empty or non-existent)
#   --cpus         threads to use
#   --attempt      current attempt number (for Nextflow retries)
#   --max-retries  maximum attempts before treating failures as soft and writing dastool.note

set -euo pipefail

assembly=""
bin_maps=""
metabat_dir=""
comebin_dir=""
semibin_dir=""
semibin_map=""
rosella_dir=""
cpus=1
attempt=0
max_retries=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --assembly)    assembly="$2";    shift 2 ;;
    --bin-maps)    bin_maps="$2";    shift 2 ;;
    --metabat-dir) metabat_dir="$2"; shift 2 ;;
    --comebin-dir) comebin_dir="$2"; shift 2 ;;
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
comebin_tsv="${tmp_dir}/comebin_bins.tsv"
rosella_tsv="${tmp_dir}/rosella_bins.tsv"
semibin_tsv="${tmp_dir}/semibin_bins.tsv"

: > "$metabat_tsv"
: > "$comebin_tsv"
: > "$rosella_tsv"
: > "$semibin_tsv"

bins_list=()

if [[ -n "$bin_maps" ]]; then
  IFS=',' read -r -a provided_maps <<< "$bin_maps"
  for map_file in "${provided_maps[@]}"; do
    [[ -s "$map_file" ]] && bins_list+=("$map_file")
  done
else
  # Legacy path kept so old calls still work while the workflow migrates to
  # normalised binner maps.
  if [[ -n "$metabat_dir" && -d "$metabat_dir" ]]; then
    for f in "$metabat_dir"/*.fa; do
      [[ -e "$f" ]] || break
      awk -v bin="${f##*/}" '
        /^>/ {
          sub(/^>/, "")
          split($0, fields, /[[:space:]]+/)
          print fields[1] "\t" bin
        }
      ' "$f" >> "$metabat_tsv"
    done
  fi

  if [[ -n "$comebin_dir" && -d "$comebin_dir" ]]; then
    for f in "$comebin_dir"/*.fa; do
      [[ -e "$f" ]] || break
      awk -v bin="${f##*/}" '
        /^>/ {
          sub(/^>/, "")
          split($0, fields, /[[:space:]]+/)
          print fields[1] "\t" bin
        }
      ' "$f" >> "$comebin_tsv"
    done
  fi

  if [[ -n "$rosella_dir" && -d "$rosella_dir" ]]; then
    for f in "$rosella_dir"/rosella_*.fna; do
      [[ -e "$f" ]] || break
      awk -v bin="${f##*/}" '
        /^>/ {
          sub(/^>/, "")
          split($0, fields, /[[:space:]]+/)
          print fields[1] "\t" bin
        }
      ' "$f" >> "$rosella_tsv"
    done
  fi

  if [[ -n "$semibin_map" && -f "$semibin_map" ]]; then
    tail -n +2 "$semibin_map" >> "$semibin_tsv"
  fi

  [[ -s "$metabat_tsv"  ]] && bins_list+=("$metabat_tsv")
  [[ -s "$comebin_tsv"  ]] && bins_list+=("$comebin_tsv")
  [[ -s "$rosella_tsv"  ]] && bins_list+=("$rosella_tsv")
  [[ -s "$semibin_tsv"  ]] && bins_list+=("$semibin_tsv")
fi

if (( ${#bins_list[@]} == 0 )); then
  msg="DASTool: no bins found for any tool"
  echo "$msg" >&2
  printf '%s\n' "$msg" > "$note_file"
  exit 0
fi

# Build comma-separated list for DAS_Tool
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

  # Special case: no high-quality bins - not a real error for our pipeline
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
