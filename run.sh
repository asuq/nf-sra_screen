#!/bin/bash
set -euo pipefail

export NXF_VER=25.04.8

# ------------------------------------------------------------------
# User settings
# ------------------------------------------------------------------
RUN_DIR='/path/to/running_dir'
NF_SRA_SCREEN='/path/to/nf-sra_screen'
ENABLE_GWDG_QOS_HELPER=false
GWDG_QOS_HELPER_OPTS=(--quiet)

QOS_HELPER_PID=""
nf_exit=0

# mamba activate /bucket/HusnikU/Conda-envs/nextflow_v25.04.8
# module load singularity

# ------------------------------------------------------------------
# Cleanup on normal exit, failure, or Ctrl+C
# ------------------------------------------------------------------
# shellcheck disable=SC2329
stop_background_process() {
  local name=$1
  local pid=$2

  if [ -n "${pid:-}" ]; then
    if kill -0 "${pid}" 2>/dev/null; then
      echo "Stopping ${name} (PID ${pid})" >&2
      kill "${pid}" 2>/dev/null || true
      wait "${pid}" 2>/dev/null || true
    fi
  fi
}

# shellcheck disable=SC2329
cleanup() {
  # Exit status of the script at the moment the trap fired
  local status=$?

  stop_background_process "gwdg_promote_2h_qos" "${QOS_HELPER_PID}"

  # Preserve the original exit status of the script
  return "${status}"
}

# Trigger cleanup on:
#  - normal exit (EXIT)
#  - Ctrl+C (INT)
#  - termination (TERM)
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

# ------------------------------------------------------------------
# Go to run directory (where .nextflow.log and execution reports will live)
# ------------------------------------------------------------------
cd "${RUN_DIR}" || exit 1

case "${ENABLE_GWDG_QOS_HELPER}" in
  true|false)
    ;;
  *)
    echo "ENABLE_GWDG_QOS_HELPER must be true or false" >&2
    exit 1
    ;;
esac

if [ "${ENABLE_GWDG_QOS_HELPER}" = true ]; then
  "${NF_SRA_SCREEN}/helpers/gwdg_promote_2h_qos.sh" \
    "${GWDG_QOS_HELPER_OPTS[@]}" \
    > gwdg_promote_2h_qos.log 2>&1 &

  QOS_HELPER_PID=$!
  echo "gwdg_promote_2h_qos pid: ${QOS_HELPER_PID}"
  printf '%s\n' "${QOS_HELPER_PID}" > gwdg_promote_2h_qos.pid
fi

# ------------------------------------------------------------------
# Run Nextflow pipeline
# ------------------------------------------------------------------
nextflow run asuq/nf-sra_screen \
  -profile <docker/singularity/local/slurm/...> \
  --sra sra.csv \
  --fastq_tsv fastq.tsv \
  --assemblers auto \
  --taxdump /path/to/ncbi_taxdump_dir \
  --uniprot_db /path/to/uniprot.dmnd \
  --taxa taxa.csv \
  --binning \
  --gtdb_ncbi_map /path/to/ncbi_vs_gtdb_xlsx_dir \
  --sandpiper_db /path/to/sandpiper_db_dir \
  --singlem_db /path/to/singlem_metapackage \
  -work-dir /lustre/path/to/nf-sra_screen_work \
  --outdir /nfs/path/to/nf-sra_screen_results \
  -resume \
|| nf_exit=$?

# This exit will trigger the EXIT trap, which runs cleanup()
exit "${nf_exit}"
