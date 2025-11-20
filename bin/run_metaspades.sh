#!/usr/bin/env bash
# Usage:
#   run_metaspades.sh \
#       --srr SRR \
#       --cpus N \
#       --memory-gb MEM \
#       --attempt A \
#       --max-retries M \
#       --reads READS.fastq.gz
#
# Produces:
#   - assembly.fasta
#   - assembly.bam
#   - assembly.bam.csi
#   - assembly.gfa
#		- spades.log
# 	- fastp.html
#   - FAIL.note (only on fatal problems)

set -euo pipefail

srr=""
attempt=0
max_retries=1
threads=1
mem_gb=4
read_files=()

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

fail() {
  local msg="$1"
  echo "$msg" >&2
  if [[ "$attempt" -lt "$max_retries" ]]; then
    exit 1
  fi
  echo "$msg" > FAIL.note
  exit 0
}


if [[ -z "$srr" ]]; then
  echo "run_metaspades.sh: --srr is required" >&2
  exit 1
fi


if (( ${#read_files[@]} < 2 )); then
  fail "metaSPAdes: paired-end reads not found"
fi

R1="${read_files[0]}"
R2="${read_files[1]}"

# fastp
if ! fastp --in1 "$R1" --in2 "$R2" --out1 "${srr}_fastp_R1.fastq.gz" --out2 "${srr}_fastp_R2.fastq.gz" \
        --thread "$threads" --length_required 50 --detect_adapter_for_pe --html fastp.html; then
  fail "metaSPAdes: fastp failed"
fi

# SPAdes
if ! spades.py -1 "${srr}_fastp_R1.fastq.gz" -2 "${srr}_fastp_R2.fastq.gz" -o '.' \
      -k 21,33,55,77,99,119,127 --meta --threads "$threads" --memory "$mem_gb"; then
  fail "metaSPAdes: assembly failed"
fi

# Rename outputs
if [[ -f scaffolds.fasta ]]; then
  mv -v scaffolds.fasta assembly.fasta
else
  mv -v contigs.fasta assembly.fasta
fi
mv -v assembly_graph_with_scaffolds.gfa assembly.gfa

# bowtie2 + samtools
if ! ( bowtie2-build -f -q --threads "$threads" assembly.fasta assembly_index \
      && bowtie2 -q --reorder --threads "$threads" --time --met-stderr --met 10 \
         -x assembly_index -1 "${srr}_fastp_R1.fastq.gz" -2 "${srr}_fastp_R2.fastq.gz" \
         | samtools sort --output-fmt BAM -@ "$threads" -o assembly.bam \
      && samtools index -c -o assembly.bam.csi -@ "$threads" assembly.bam ); then
  fail "metaSPAdes: mapping/indexing failed"
fi
