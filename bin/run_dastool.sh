#!/usr/bin/env bash
# Run DAS_Tool to integrate bins from MetaBAT, CONCOCT, SemiBin, and Rosella

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

if [[ -z "$assembly" || -z "$metabat_dir" || -z "$concoct_dir" || -z "$semibin_dir" || -z "$semibin_map" || -z "$rosella_dir" ]]; then
  echo "run_dastool.sh: missing required arguments" >&2
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

tmp_dir="tmp_dastool"
final_dir="dastool"

rm -rf "$tmp_dir" "$final_dir"
mkdir -p "$tmp_dir"

shopt -s nullglob

# Build contig2bin TSVs
metabat_tsv="${tmp_dir}/metabat_bins.tsv"
concoct_tsv="${tmp_dir}/concoct_bins.tsv"
rosella_tsv="${tmp_dir}/rosella_bins.tsv"
semibin_tsv="${tmp_dir}/semibin_bins.tsv"

: > "$metabat_tsv"
: > "$concoct_tsv"
: > "$rosella_tsv"

for f in "$metabat_dir"/*.fa; do
  [[ -e "$f" ]] || break
  grep "^>" "$f" | cut -f1 | sed 's/>//' | sed "s/$/\t${f##*/} /" >> "$metabat_tsv"
done

for f in "$concoct_dir"/*.fa; do
  [[ -e "$f" ]] || break
  grep "^>" "$f" | sed 's/>.*/&\t'${f##*/}' /' | sed 's/>//' >> "$concoct_tsv"
done

for f in "$rosella_dir"/rosella_*.fna; do
  [[ -e "$f" ]] || break
  grep "^>" "$f" | sed 's/>.*/&\t'${f##*/}' /' | sed 's/>//' >> "$rosella_tsv"
done

if [[ ! -f "$semibin_map" ]]; then
  fail "DASTool: SemiBin contig_bins.tsv not found at $semibin_map"
fi

tail -n +2 "$semibin_map" > "$semibin_tsv"

if [[ ! -s "$metabat_tsv" && ! -s "$concoct_tsv" && ! -s "$rosella_tsv" && ! -s "$semibin_tsv" ]]; then
  fail "DASTool: no bins found for any tool"
fi

if ! DAS_Tool \
      --bins "${metabat_tsv},${concoct_tsv},${rosella_tsv},${semibin_tsv}" \
      --contigs "$assembly" \
      --outputbasename "${tmp_dir}/dastool" \
      --write_bin_evals \
      --write_bins \
      --threads "$cpus"; then
  fail "DASTool: refinement failed"
fi

if [[ ! -d "${tmp_dir}/dastool_DASTool_bins" ]]; then
  fail "DASTool: dastool_DASTool_bins directory not found"
fi

mv "${tmp_dir}/dastool_DASTool_bins" "$final_dir"
