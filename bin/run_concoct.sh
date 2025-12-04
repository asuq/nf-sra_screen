#!/usr/bin/env bash
# Run CONCOCT binning on an assembly + BAM
#
# Args:
#   --assembly     assembly fasta
#   --bam          BAM mapped against the assembly
#   --cpus         threads to use
#   --attempt      current attempt number (for Nextflow retries)
#   --max-retries  maximum attempts before treating failures as soft and writing concoct.note

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
      echo "run_concoct.sh: unknown argument: $1" >&2
      exit 1 ;;
  esac
done

if [[ -z "$assembly" || -z "$bam" ]]; then
  echo "run_concoct.sh: missing --assembly or --bam" >&2
  exit 1
fi

tmp_dir="tmp_concoct"
final_dir="concoct"
note_file="concoct.note"

rm -rf "$tmp_dir" "$final_dir" "$note_file"
mkdir -p "$tmp_dir"
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

# Cut contigs into ~10 kb chunks
if ! cut_up_fasta.py "$assembly" -c 10000 -o 0 --merge_last \
     -b "${tmp_dir}/contigs_10k.bed" > "${tmp_dir}/contigs_10k.fa"; then
  fail "Concoct: cut_up_fasta.py failed"
fi

# Coverage table from BAM
if ! concoct_coverage_table.py "${tmp_dir}/contigs_10k.bed" "$bam" > "${tmp_dir}/coverage_table.tsv"; then
  fail "Concoct: concoct_coverage_table.py failed"
fi

# Run CONCOCT
if ! concoct --composition_file "${tmp_dir}/contigs_10k.fa" \
             --coverage_file "${tmp_dir}/coverage_table.tsv" \
             --basename "${tmp_dir}/concoct" \
             --threads "$cpus"; then
  fail "Concoct: concoct failed"
fi

# Merge clusters back to full contigs
if ! merge_cutup_clustering.py "${tmp_dir}/concoct_clustering_gt1000.csv" \
     > "${tmp_dir}/concoct_clustering_merged.csv"; then
  fail "Concoct: merge_cutup_clustering.py failed"
fi


# Extract bins
if ! extract_fasta_bins.py "$assembly" "${tmp_dir}/concoct_clustering_merged.csv" \
       --output_path "$final_dir"; then
  fail "Concoct: extract_fasta_bins.py failed"
fi
