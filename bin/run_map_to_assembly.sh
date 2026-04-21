#!/usr/bin/env bash
# Usage:
#   run_map_to_assembly.sh \
#     --read-type READ_TYPE \
#     --assembly ASSEMBLY_FASTA \
#     --cpus N \
#     --attempt A \
#     --max-retries M \
#     --reads READS.fastq.gz
#
# Produces:
#   - assembly.bam
#   - assembly.bam.csi
#   - FAIL.note (only on fatal problems)

set -euo pipefail

read_type=""
assembly=""
cpus=1
attempt=0
max_retries=1
read_files=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --reads)
      shift
      while [[ $# -gt 0 ]]; do
        case "$1" in
          --read-type|--assembly|--cpus|--attempt|--max-retries)
            break
            ;;
          *)
            read_files+=("$1")
            shift
            ;;
        esac
      done
      ;;
    --read-type)   read_type="$2"; shift 2 ;;
    --assembly)    assembly="$2"; shift 2 ;;
    --cpus)        cpus="$2"; shift 2 ;;
    --attempt)     attempt="$2"; shift 2 ;;
    --max-retries) max_retries="$2"; shift 2 ;;
    *)
      echo "run_map_to_assembly.sh: unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

if [[ -z "$read_type" || -z "$assembly" ]]; then
  echo "run_map_to_assembly.sh: missing --read-type or --assembly" >&2
  exit 1
fi

if [[ ${#read_files[@]} -eq 0 ]]; then
  echo "run_map_to_assembly.sh: missing --reads" >&2
  exit 1
fi

rm -f FAIL.note assembly.bam assembly.bam.csi

fail() {
  local msg="$1"
  echo "$msg" >&2
  if [[ "$attempt" -lt "$max_retries" ]]; then
    exit 1
  fi
  printf '%s\n' "$msg" > FAIL.note
  exit 0
}

map_short_reads() {
  if (( ${#read_files[@]} == 1 )); then
    bowtie2-build -f -q --threads "$cpus" "$assembly" assembly_index
    bowtie2 -q --reorder --threads "$cpus" --time --met-stderr --met 10 \
      -x assembly_index \
      -U "${read_files[0]}" \
      | samtools sort --output-fmt BAM -@ "$cpus" -o assembly.bam
    return
  fi

  if (( ${#read_files[@]} == 2 )); then
    bowtie2-build -f -q --threads "$cpus" "$assembly" assembly_index
    bowtie2 -q --reorder --threads "$cpus" --time --met-stderr --met 10 \
      -x assembly_index \
      -1 "${read_files[0]}" \
      -2 "${read_files[1]}" \
      | samtools sort --output-fmt BAM -@ "$cpus" -o assembly.bam
    return
  fi

  fail "Mapping (short): expected one or two FASTQ files"
}

map_long_reads() {
  local preset="$1"
  minimap2 -ax "$preset" -I 20G -t "$cpus" "$assembly" "${read_files[@]}" \
    | samtools sort --output-fmt BAM -@ "$cpus" -o assembly.bam
}

case "${read_type,,}" in
  short)
    if ! map_short_reads; then
      fail "Mapping (short): alignment failed"
    fi
    ;;
  nanopore)
    if ! map_long_reads "map-ont"; then
      fail "Mapping (nanopore): alignment failed"
    fi
    ;;
  pacbio)
    if ! map_long_reads "map-pb"; then
      fail "Mapping (pacbio): alignment failed"
    fi
    ;;
  hifi)
    if ! map_long_reads "map-hifi"; then
      fail "Mapping (hifi): alignment failed"
    fi
    ;;
  *)
    echo "run_map_to_assembly.sh: unsupported read type '${read_type}'" >&2
    exit 1
    ;;
esac

if ! samtools index -c -o assembly.bam.csi -@ "$cpus" assembly.bam; then
  fail "Mapping (${read_type}): indexing failed"
fi
