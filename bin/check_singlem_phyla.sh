#!/usr/bin/env bash
# check_singlem_phyla.sh
# $1: path to singlem_taxonomic_profile.tsv
# $2: path to phyla_to_check.txt
# $3: output path for singlem_output.tsv
# $4: task.attempt (provided by Nextflow)
# $5: params.max_retries (provided by Nextflow)

function fail() {
  local msg="$1"
  echo "$msg" >&2
  if [[ $4 -lt $5 ]]; then
    exit 1
  fi
  echo "$msg" > FAIL.note
  exit 0
}


rc=0
if check_singlem_phyla.py -i "$1" -p "$2" -o "$3"; then
  shopt -s nullglob
  mkdir -p reads.ok
  for f in *.f*q*; do
    ln -sf "../$f" "reads.ok/$(basename "$f")"
  done
  shopt -u nullglob
  rc=0

else
  rc=$?
fi

case "$rc" in
  0)
    : # all good; emit reads.ok/*
    ;;
  1)
    fail "SingleM phylum check internal error"
    ;;
  2)
    fail "No target phyla detected by SingleM"
    ;;
esac
