#!/bin/bash
set -euo pipefail

export NXF_VER=25.04.8

# ------------------------------------------------------------------
# User settings
# ------------------------------------------------------------------
RUN_DIR='/path/to/running_dir'
DEST_DIR='/path/to/destination_dir'
INTERVAL_MIN=10
NF_SRA_SCREEN='/path/to/nf-sra_screen'

WATCH_PID=""

# mamba activate /bucket/HusnikU/Conda-envs/nextflow_v25.04.8
# module load singularity

# ------------------------------------------------------------------
# Cleanup on normal exit, failure, or Ctrl+C
# ------------------------------------------------------------------
cleanup() {
  # Exit status of the script at the moment the trap fired
  local status=$?

  if [ -n "${WATCH_PID:-}" ]; then
    if kill -0 "${WATCH_PID}" 2>/dev/null; then
      echo "Stopping watch_and_transfer (PID ${WATCH_PID})" >&2
      # Try to terminate gracefully; ignore errors if it is already gone
      kill "${WATCH_PID}" 2>/dev/null || true
      # Reap it so it does not become a zombie
      wait "${WATCH_PID}" 2>/dev/null || true
    fi
  fi

  # Preserve the original exit status of the script
  return "${status}"
}

# Trigger cleanup on:
#  - normal exit (EXIT)
#  - Ctrl+C (INT)
#  - termination (TERM)
trap cleanup EXIT INT TERM

# ------------------------------------------------------------------
# Go to run directory (where .nextflow.log, work/, output/ will live)
# ------------------------------------------------------------------
cd "${RUN_DIR}" || exit 1

# ------------------------------------------------------------------
# Start watcher in background
# ------------------------------------------------------------------
mkdir -p "${DEST_DIR}"

"${NF_SRA_SCREEN}/bin/watch_and_transfer.sh" \
  "${RUN_DIR}" \
  "${DEST_DIR}" \
  "${INTERVAL_MIN}" \
  > watch_and_transfer.log 2>&1 &

WATCH_PID=$!
echo "watch_and_transfer pid: ${WATCH_PID}"
printf '%s\n' "${WATCH_PID}" > watch_and_transfer.pid

# ------------------------------------------------------------------
# Run Nextflow pipeline (do not let set -e hide its exit code)
# ------------------------------------------------------------------
set +e
nextflow run asuq/nf-sra_screen \
  -profile <docker/singularity/local/slurm/...> \
  --sra sra.csv \
  --fastq_tsv fastq.tsv \
  --taxdump /path/to/ncbi_taxdump_dir \
  --uniprot_db /path/to/uniprot.dmnd \
  --taxa taxa.csv \
  --binning \
  --gtdb_ncbi_map /path/to/ncbi_vs_gtdb_xlsx_dir \
  --sandpiper_db /path/to/sandpiper_db_dir \
  --singlem_db /path/to/singlem_metapackage \
  --outdir nf-sra_screen_results \
  -resume
nf_exit=$?
set -e

# This exit will trigger the EXIT trap, which runs cleanup()
exit "${nf_exit}"
