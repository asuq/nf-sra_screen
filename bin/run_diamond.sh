#!/usr/bin/env bash
# Usage:
#   run_diamond.sh \
#       --assembly assembly.fasta \
#       --db uniprot.dmnd \
#       --cpus N \
#       --attempt A \
#       --max-retries M
#
# Produces:
#   - assembly_vs_uniprot.tsv
#   - FAIL.note (only on fatal problems)

set -euo pipefail

assembly=""
db=""
attempt=0
max_retries=1
cpus=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --assembly)    assembly="$2"; shift 2 ;;
    --db)          db="$2"; shift 2 ;;
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

if ! diamond blastx --sensitive --query "$assembly" \
      --out "assembly_vs_uniprot.tsv" --db "$db" \
      --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
      --verbose --threads "$cpus" --evalue 1e-25 --max-target-seqs 5; then
  fail "Diamond: run failed"
fi
