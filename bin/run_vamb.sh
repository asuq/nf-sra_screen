#!/usr/bin/env bash
# Run VAMB binning on an assembly + BAM.
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
      echo "run_vamb.sh: unknown argument: $1" >&2
      exit 1 ;;
  esac
done

if [[ -z "$assembly" || -z "$bam" ]]; then
  echo "run_vamb.sh: missing --assembly or --bam" >&2
  exit 1
fi

tmp_dir="tmp_vamb"
bam_dir="vamb_bam"
final_dir="vamb"
note_file="vamb.note"
map_file="vamb.contig2bin.tsv"

rm -rf "$tmp_dir" "$bam_dir" "$final_dir" "$note_file" "$map_file"
mkdir -p "$bam_dir"
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

ln -s "$(realpath "$bam")" "${bam_dir}/$(basename "$bam")"
if [[ -e "${bam}.bai" ]]; then
  ln -s "$(realpath "${bam}.bai")" "${bam_dir}/$(basename "${bam}.bai")"
fi
if [[ -e "${bam}.csi" ]]; then
  ln -s "$(realpath "${bam}.csi")" "${bam_dir}/$(basename "${bam}.csi")"
fi

if ! vamb bin default \
      --outdir "$tmp_dir" \
      --fasta "$assembly" \
      --bamdir "$bam_dir" \
      --minfasta 1 \
      -p "$cpus"; then
  fail "VAMB: bin default failed"
fi

cluster_file=""
if [[ -f "${tmp_dir}/vae_clusters_split.tsv" ]]; then
  cluster_file="${tmp_dir}/vae_clusters_split.tsv"
elif [[ -f "${tmp_dir}/vae_clusters_unsplit.tsv" ]]; then
  cluster_file="${tmp_dir}/vae_clusters_unsplit.tsv"
else
  fail "VAMB: cluster TSV not found"
fi

awk 'BEGIN { FS = OFS = "\t" }
  tolower($1) == "clustername" && tolower($2) == "contigname" {
    next
  }
  NF >= 2 && $1 != "" && $2 != "" {
    print $2, $1
  }
' "$cluster_file" > "$map_file"

mkdir -p "$final_dir"
if [[ -d "${tmp_dir}/bins" ]]; then
  shopt -s nullglob
  for f in "${tmp_dir}/bins"/*; do
    [[ -f "$f" ]] || continue
    cp "$f" "$final_dir"/
  done
  shopt -u nullglob
fi
