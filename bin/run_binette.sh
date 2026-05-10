#!/usr/bin/env bash
# Run Binette to refine normalised binner maps.
#
# Args:
#   --assembly     assembly fasta
#   --bin-maps     comma-separated contig-to-bin TSVs from selected binners
#   --checkm2-db   CheckM2 DIAMOND database
#   --cpus         threads to use
#   --attempt      current attempt number (for Nextflow retries)
#   --max-retries  maximum attempts before treating failures as soft

set -euo pipefail

assembly=""
bin_maps=""
checkm2_db=""
cpus=1
attempt=0
max_retries=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --assembly)    assembly="$2";    shift 2 ;;
    --bin-maps)    bin_maps="$2";    shift 2 ;;
    --checkm2-db)  checkm2_db="$2";  shift 2 ;;
    --cpus)        cpus="$2";        shift 2 ;;
    --attempt)     attempt="$2";     shift 2 ;;
    --max-retries) max_retries="$2"; shift 2 ;;
    *)
      echo "run_binette.sh: unknown argument: $1" >&2
      exit 1 ;;
  esac
done

if [[ -z "$assembly" || -z "$checkm2_db" ]]; then
  echo "run_binette.sh: --assembly and --checkm2-db are required" >&2
  exit 1
fi

resolve_checkm2_db() {
  # Binette forwards this path to DIAMOND, which needs the database file.
  local db_path="$1"
  local default_db="${db_path}/uniref100.KO.1.dmnd"

  if [[ -d "$db_path" ]]; then
    if [[ -f "$default_db" ]]; then
      printf '%s\n' "$default_db"
      return 0
    fi
    echo "run_binette.sh: CheckM2 DB directory is missing uniref100.KO.1.dmnd: ${db_path}" >&2
    exit 1
  fi

  printf '%s\n' "$db_path"
}

resolved_checkm2_db="$(resolve_checkm2_db "$checkm2_db")"

tmp_dir="tmp_binette"
final_dir="binette"
note_file="binette.note"

rm -rf "$tmp_dir" "$final_dir" "$note_file"
mkdir -p "$tmp_dir"
: > "$note_file"

record_soft_failure() {
  local msg="$1"
  rm -rf "$final_dir"
  mkdir -p "$final_dir"
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

  record_soft_failure "Binette: unexpected failure (exit ${exit_code})"
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

bins_list=()
if [[ -n "$bin_maps" ]]; then
  IFS=',' read -r -a provided_maps <<< "$bin_maps"
  for map_file in "${provided_maps[@]}"; do
    [[ -s "$map_file" ]] && bins_list+=("$map_file")
  done
fi

if (( ${#bins_list[@]} == 0 )); then
  msg="Binette: no non-empty binner maps found"
  echo "$msg" >&2
  printf '%s\n' "$msg" > "$note_file"
  mkdir -p "$final_dir"
  exit 0
fi

if binette \
      --contig2bin_tables "${bins_list[@]}" \
      --contigs "$assembly" \
      --checkm2_db "$resolved_checkm2_db" \
      --outdir "$tmp_dir" \
      --prefix binette \
      --threads "$cpus"; then
  :
else
  fail "Binette: refinement failed" "$?"
fi

if [[ ! -d "${tmp_dir}/final_bins" && ! -f "${tmp_dir}/final_contig_to_bin.tsv" ]]; then
  fail "Binette: final output not found"
fi

mv "$tmp_dir" "$final_dir"
