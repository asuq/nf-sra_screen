#!/usr/bin/env bash
# Usage:
#   run_append_summary.sh \
#       --outdir OUTDIR \
#       --sra SRA \
#       --srr SRR \
#       --platform PLATFORM \
#       --model MODEL \
#       --strategy STRATEGY \
#       --assembler ASSEMBLER \
#       --summary-csv summary.csv \
#       --note NOTE

set -euo pipefail

outdir=""
sra=""
srr=""
platform=""
model=""
strategy=""
assembler=""
summary_csv=""
note=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir)      outdir="$2"; shift 2 ;;
    --sra)         sra="$2"; shift 2 ;;
    --srr)         srr="$2"; shift 2 ;;
    --platform)    platform="$2"; shift 2 ;;
    --model)       model="$2"; shift 2 ;;
    --strategy)    strategy="$2"; shift 2 ;;
    --assembler)   assembler="$2"; shift 2 ;;
    --summary-csv) summary_csv="$2"; shift 2 ;;
    --note)        note="$2"; shift 2 ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1 ;;
  esac
done

if [[ -z "${outdir}" || -z "${sra}" || -z "${summary_csv}" ]]; then
  echo "append_summary.sh: missing required parameters" >&2
  exit 1
fi

OUT_TSV="${outdir}/summary.tsv"
mkdir -p "${outdir}"

# counts = comma-joined n_identifiers in the order of rows in summary.csv
COUNTS="$(
  awk -F',' 'NR>1{print $3}' "${summary_csv}" | paste -sd, -
)"

# note left blank for successful path
LINE="${sra}\t${srr}\t${platform}\t${model}\t${strategy}\t${assembler}\t${COUNTS}\t${note}"

{
  flock 200
  if [[ ! -s "${OUT_TSV}" ]]; then
    printf 'sra\tsrr\tplatform\tmodel\tstrategy\tassembler\tcounts\tnote\n' > "${OUT_TSV}"
  fi
  printf '%b\n' "${LINE}" >> "${OUT_TSV}"
} 200> "${OUT_TSV}.lock"

ln -sf "${OUT_TSV}" .
