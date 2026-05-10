#!/usr/bin/env bash
# Run SemiBin2 single_easy_bin on an assembly + BAM
# Uses DIAMONDDB from --diamond-db (same DB as DIAMOND process)
#
# Args:
#   --assembly     assembly fasta
#   --bam          BAM mapped against the assembly
#   --diamond-db   DIAMOND reference database
#   --read-type    read type: short, nanopore, pacbio, or hifi
#   --environment  SemiBin2 pretrained environment
#   --cpus         threads to use
#   --attempt      current attempt number (for Nextflow retries)
#   --max-retries  maximum attempts before treating failures as soft and writing semibin.note

set -euo pipefail

assembly=""
bam=""
diamond_db=""
read_type=""
environment="global"
cpus=1
attempt=0
max_retries=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --assembly) assembly="$2"; shift 2 ;;
    --bam) bam="$2"; shift 2 ;;
    --diamond-db) diamond_db="$2"; shift 2 ;;
    --read-type) read_type="$2"; shift 2 ;;
    --environment) environment="$2"; shift 2 ;;
    --cpus) cpus="$2"; shift 2 ;;
    --attempt) attempt="$2"; shift 2 ;;
    --max-retries) max_retries="$2"; shift 2 ;;
    *)
      echo "run_semibin.sh: unknown argument: $1" >&2
      exit 1 ;;
  esac
done

if [[ -z "$assembly" || -z "$bam" || -z "$diamond_db" || -z "$read_type" ]]; then
  echo "run_semibin.sh: missing --assembly, --bam, --diamond-db, or --read-type" >&2
  exit 1
fi

export DIAMONDDB="$diamond_db"

read_type_lc="$(printf '%s' "$read_type" | tr '[:upper:]' '[:lower:]')"
environment_lc="$(printf '%s' "$environment" | tr '[:upper:]' '[:lower:]')"

case "$read_type_lc" in
  short)
    sequencing_type="short_read"
    ;;
  nanopore|pacbio|hifi)
    sequencing_type="long_read"
    ;;
  *)
    echo "run_semibin.sh: unsupported read type '${read_type}'" >&2
    exit 1
    ;;
esac

case "$environment_lc" in
  human_gut|dog_gut|ocean|soil|cat_gut|human_oral|mouse_gut|pig_gut|built_environment|wastewater|chicken_caecum|global)
    ;;
  *)
    echo "run_semibin.sh: unsupported SemiBin environment '${environment}'" >&2
    exit 1
    ;;
esac

tmp_dir="tmp_semibin"
final_dir="semibin"
note_file="semibin.note"
map_file="semibin.contig2bin.tsv"

rm -rf "$tmp_dir" "$final_dir" "$note_file" "$map_file"
mkdir -p "$tmp_dir" "$final_dir"
: > "$note_file"
: > "$map_file"

record_soft_failure() {
  local msg="$1"
  rm -rf "$final_dir"
  mkdir -p "$final_dir"
  touch "$final_dir/contig_bins.tsv"
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

  record_soft_failure "SemiBin2: unexpected failure (exit ${exit_code})"
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

if SemiBin2 single_easy_bin \
      --input-fasta "$assembly" \
      --input-bam "$bam" \
      --environment "$environment_lc" \
      --sequencing-type "$sequencing_type" \
      --engine cpu \
      --output "$tmp_dir" \
      --threads "$cpus"; then
  :
else
  fail "SemiBin2: single_easy_bin failed" "$?"
fi

if [[ ! -f "${tmp_dir}/contig_bins.tsv" ]]; then
  fail "SemiBin2: contig_bins.tsv not found"
fi

cp "${tmp_dir}/contig_bins.tsv" "${final_dir}/contig_bins.tsv"
tail -n +2 "${tmp_dir}/contig_bins.tsv" > "$map_file"

shopt -s nullglob

# Collect any bins produced (gzipped or not)
bins=( "${tmp_dir}/output_bins"/*.gz "${tmp_dir}/output_bins"/* )

if [[ "${#bins[@]}" -ne 0 ]]; then
  # Copy bins to final directory
  for f in "${bins[@]}"; do
    cp "$f" "$final_dir"/
  done

  # Parallel gunzip, if needed
  gz_files=( "$final_dir"/*.gz )
  if [[ "${#gz_files[@]}" -gt 0 ]]; then
    printf '%s\0' "${gz_files[@]}" \
      | xargs -0 -n 1 -P "${cpus}" gunzip -v \
      || fail "SemiBin2: gunzip failed" "$?"
  fi
fi

shopt -u nullglob
