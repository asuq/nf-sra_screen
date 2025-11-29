#!/usr/bin/env bash
# Usage:
#   run_blobtools.sh \
#       --assembly assembly.fasta \
#       --bam assembly.bam \
#       --csi assembly.bam.csi \
#       --blast assembly_vs_uniprot.tsv \
#       --taxdump /path/to/taxdump \
#       --cpus N \
#       --attempt A \
#       --max-retries M
# Prod
#   - blobtools/ (directory with blobtools outputs)
#   - blobtools.csv
#   - FAIL.note (only on fatal problems)

set -euo pipefail

assembly=""
bam=""
csi=""
blast=""
taxdump=""
cpus=1
attempt=0
max_retries=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --assembly)    assembly="$2"; shift 2 ;;
    --bam)         bam="$2"; shift 2 ;;
    --csi)         csi="$2"; shift 2 ;;
    --blast)       blast="$2"; shift 2 ;;
    --taxdump)     taxdump="$2"; shift 2 ;;
    --cpus)        cpus="$2"; shift 2 ;;
    --attempt)     attempt="$2"; shift 2 ;;
    --max-retries) max_retries="$2"; shift 2 ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

fail() {
  local msg="$1"
  echo "$msg" >&2
  if [[ "$attempt" -lt "$max_retries" ]]; then
    exit 1
  fi
  echo "$msg" > FAIL.note
  exit 0
}

if ! blobtools create --fasta "$assembly" --cov "$bam" \
      --hits "$blast" --taxdump "$taxdump" --threads "$cpus" 'blobtools'; then
  fail "Blobtools: failed to create blobdir"
fi

if ! blobtools filter --table 'blobtools.tsv' \
      --table-fields gc,length,ncount,assembly_cov,bestsumorder_superkingdom,bestsumorder_kingdom,bestsumorder_phylum,bestsumorder_class,bestsumorder_order,bestsumorder_family,bestsumorder_genus,bestsumorder_species 'blobtools'; then
  fail "Blobtools: failed to filter result"
fi

if ! blobtools view --format png --out 'blobtools' --plot 'blobtools'; then
	fail "Blobtools: failed to generate plots"
fi

# TSV → CSV
header=$(head -n 1 'blobtools.tsv' | cut -f2- | sed 's/bestsumorder_//g' | sed "s/assembly_cov/coverage/" | tr '\t' ',')
data=$(tail -n +2 'blobtools.tsv' | cut -f2- | tr '\t' ',')
{ echo "${header}"; echo "${data}"; } > 'blobtools.csv'
