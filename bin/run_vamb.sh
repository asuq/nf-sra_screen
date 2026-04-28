#!/usr/bin/env bash
# Run VAMB binning on an assembly + BAM.
#
# Args:
#   --assembly    assembly fasta
#   --bam         BAM mapped against the assembly
#   --cpus        threads to use
#   --attempt     current attempt number (for Nextflow retries)
#   --max-retries maximum attempts before treating failures as soft
#   --cuda        pass VAMB's CUDA flag for model training and clustering
#   --require-cuda fail before binning if PyTorch cannot see CUDA

set -euo pipefail

assembly=""
bam=""
cpus=1
attempt=0
max_retries=1
use_cuda=false
require_cuda=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --assembly) assembly="$2"; shift 2 ;;
    --bam) bam="$2"; shift 2 ;;
    --cpus) cpus="$2"; shift 2 ;;
    --attempt) attempt="$2"; shift 2 ;;
    --max-retries) max_retries="$2"; shift 2 ;;
    --cuda) use_cuda=true; shift ;;
    --require-cuda) require_cuda=true; shift ;;
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

require_pytorch_cuda() {
  python - <<'PY'
"""Fail clearly unless PyTorch can access a CUDA device."""
import importlib.util
import sys


def main() -> int:
    """Check PyTorch CUDA visibility."""
    if importlib.util.find_spec("torch") is None:
        print("VAMB GPU mode requires PyTorch, but torch is not installed", file=sys.stderr)
        return 1

    import torch

    if not torch.cuda.is_available():
        print("VAMB GPU mode requested, but CUDA is not visible to PyTorch", file=sys.stderr)
        return 1

    print(
        f"VAMB GPU preflight: CUDA visible to PyTorch "
        f"(cuda={torch.version.cuda}, devices={torch.cuda.device_count()})",
        file=sys.stderr,
    )
    return 0


raise SystemExit(main())
PY
}

if [[ "$require_cuda" == true ]]; then
  require_pytorch_cuda || fail "VAMB: GPU mode requested but CUDA preflight failed"
fi

ln -s "$(realpath "$bam")" "${bam_dir}/$(basename "$bam")"
if [[ -e "${bam}.bai" ]]; then
  ln -s "$(realpath "${bam}.bai")" "${bam_dir}/$(basename "${bam}.bai")"
fi
if [[ -e "${bam}.csi" ]]; then
  ln -s "$(realpath "${bam}.csi")" "${bam_dir}/$(basename "${bam}.csi")"
fi

vamb_args=(
  vamb bin default
  --outdir "$tmp_dir"
  --fasta "$assembly"
  --bamdir "$bam_dir"
  --minfasta 1
  -p "$cpus"
)
if [[ "$use_cuda" == true ]]; then
  vamb_args+=(--cuda)
fi

if ! "${vamb_args[@]}"; then
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
