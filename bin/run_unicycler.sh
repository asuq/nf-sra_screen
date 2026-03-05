#!/usr/bin/env bash
# Usage:
#   run_unicycler.sh \
#     --srr SRR \
#     --cpus N \
#     --memory-gb MEM \
#     --attempt A \
#     --max-retries M \
#     --reads READS.fastq.gz
#
# Produces:
#   - assembly.fasta
#   - assembly.bam
#   - assembly.bam.csi
#   - assembly.gfa
#   - spades.log
#   - fastp.html
#   - FAIL.note (only on fatal problems)

set -euo pipefail

srr=""
attempt=0
max_retries=1
threads=1
mem_gb=4
read_files=()


cleanup() {
  local exit_status=$?

  # Never let cleanup failures replace the real script status.
  set +e

  case "${work_dir:-}" in
    ""|"/")
      printf 'Refusing cleanup in unsafe directory: %q\n' "${work_dir:-<unset>}" >&2
      return "$exit_status"
      ;;
  esac

  find "$work_dir" \
    -mindepth 1 \
    -maxdepth 1 \
    -type d \
    -exec rm -rf -- {} +

  return "$exit_status"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --reads)
      shift
      while [[ $# -gt 0 ]]; do
        # Stop when we hit the next flag
        case "$1" in
          --sandpiper-decision|--valid-taxa|--singlem-db|--cpus|--attempt|--max-retries)
            break
            ;;
          *)
            read_files+=("$1")
            shift
            ;;
        esac
      done
      ;;
    --srr)         srr="$2"; shift 2 ;;
    --cpus)        threads="$2"; shift 2 ;;
    --memory-gb)   mem_gb="$2"; shift 2 ;;
    --attempt)     attempt="$2"; shift 2 ;;
    --max-retries) max_retries="$2"; shift 2 ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

if [[ -z "$srr" ]]; then
  echo "run_unicycler.sh: --srr is required" >&2
  exit 1
fi


if (( ${#read_files[@]} < 2 )); then
  fail "Unicycler: paired-end reads not found"
fi

readonly work_dir=$PWD
trap cleanup EXIT

fail() {
  local msg="$1"
  echo "$msg" >&2
  if [[ "$attempt" -lt "$max_retries" ]]; then
    exit 1
  fi
  echo "$msg" > FAIL.note
  exit 0
}

R1="${read_files[0]}"
R2="${read_files[1]}"

# fastp
if ! fastp --in1 "${R1}" --in2 "${R2}" --out1 "${srr}_fastp_R1.fastq.gz" --out2 "${srr}_fastp_R2.fastq.gz" \
        --thread "${threads}" --length_required 50 --detect_adapter_for_pe --qualified_quality_phred 30 --html fastp.html; then
  fail "Unicycler: fastp failed"
fi

# Unicycler
if ! unicycler -1 "${srr}_fastp_R1.fastq.gz" -2 "${srr}_fastp_R2.fastq.gz" -o '.' \
      --threads "${threads}" --verbosity 2 ; then
  fail "Unicycler: assembly failed"
fi

# bowtie2 + samtools
if ! ( bowtie2-build -f -q --threads "$threads" assembly.fasta assembly_index \
      && bowtie2 -q --reorder --threads "$threads" --time --met-stderr --met 10 \
         -x assembly_index -1 "${srr}_fastp_R1.fastq.gz" -2 "${srr}_fastp_R2.fastq.gz" \
         | samtools sort --output-fmt BAM -@ "$threads" -o assembly.bam \
      && samtools index -c -o assembly.bam.csi -@ "$threads" assembly.bam ); then
  fail "Unicycler: mapping/indexing failed"
fi
