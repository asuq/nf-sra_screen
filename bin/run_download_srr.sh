#!/usr/bin/env bash
# Usage:
#   run_download_srr.sh \
#     --srr SRR \
#     --platform PLATFORM \
#     --assembler ASM \
#     --cpus N \
#     --attempt A \
#     --max-retries M
#
# Produces:
#   - *.fastq.gz
#   - assembler.txt
#   - FAIL.note (only on fatal problems)

set -euo pipefail

srr=""
platform=""
assembler=""
cpus=1
attempt=0
max_retries=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --srr) srr="$2"; shift 2 ;;
    --platform) platform="$2"; shift 2 ;;
    --assembler) assembler="$2"; shift 2 ;;
    --cpus) cpus="$2"; shift 2 ;;
    --attempt) attempt="$2"; shift 2 ;;
    --max-retries) max_retries="$2"; shift 2 ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1 ;;
  esac
done

if [[ -z "${srr}" || -z "${platform}" ]]; then
  echo "run_download_srr.sh: missing required arguments (--srr, --platform)" >&2
  exit 1
fi

# Internal retry count (inside a single Nextflow task attempt)
INTERNAL_RETRIES=5

# Aria2 params defaults
ARIA2_J_CAP="${ARIA2_J_CAP:-2}"          # max concurrent files per SRR
ARIA2_SPLIT="${ARIA2_SPLIT:-4}"          # segments per file
ARIA2_MAX_CONN="${ARIA2_MAX_CONN:-4}"    # connections per server per file
ARIA2_TIMEOUT="${ARIA2_TIMEOUT:-60}"
ARIA2_CONNECT_TIMEOUT="${ARIA2_CONNECT_TIMEOUT:-30}"
ARIA2_RETRY_WAIT="${ARIA2_RETRY_WAIT:-5}"
ARIA2_MAX_TRIES="${ARIA2_MAX_TRIES:-2}"

log() {
  printf '[%s] %s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" "$*" >&2
}

backoff_sleep() {
  # Exponential backoff with a small deterministic jitter
  # try_idx starts at 1
  local try_idx="$1"
  local base=$(( 2 ** (try_idx - 1) ))
  local jitter=$(( (try_idx * 7) % 5 ))
  local sleep_s=$(( base + jitter ))
  if (( sleep_s > 60 )); then sleep_s=60; fi
  sleep "${sleep_s}"
}

cleanup_partial() {
  rm -rf tmp_srr 2>/dev/null || true
  rm -f ./*.fastq ./*.fq ./*.fasta ./*.fa 2>/dev/null || true
  rm -f ./*.aria2 2>/dev/null || true
}

normalise_ena_url() {
  # ENA often returns paths like "ftp.sra.ebi.ac.uk/vol1/fastq/..."
  local u="$1"
  if [[ -z "${u}" ]]; then
    printf '%s' ""
    return 0
  fi
  if [[ "${u}" == http://* || "${u}" == https://* || "${u}" == ftp://* ]]; then
    printf '%s' "${u}"
    return 0
  fi
  printf 'https://%s' "${u}"
}

is_fastq_name() {
  local f="$1"
  case "${f}" in
    *.fastq|*.fq|*.fastq.gz|*.fq.gz|*.fastq.bz2|*.fq.bz2) return 0 ;;
    *) return 1 ;;
  esac
}

md5_ok() {
  local file="$1"
  local expected="$2"
  if [[ -z "${expected}" ]]; then
    return 1
  fi
  if [[ ! -s "${file}" ]]; then
    return 1
  fi
  local got
  got="$(md5sum "${file}" | awk '{print $1}')"
  [[ "${got}" == "${expected}" ]]
}

aria2_download_with_md5() {
  # Args: url out_file expected_md5
  local url="$1"
  local out_file="$2"
  local expected_md5="$3"

  local tmp_input
  tmp_input="$(mktemp aria2_input.XXXXXX.txt)"

  # aria2 input file format: URI line, then indented option lines apply to it.
  # We also request checksum checking; still do an explicit md5sum afterwards.
  {
    printf '%s\n' "${url}"
    printf '  out=%s\n' "${out_file}"
    if [[ -n "${expected_md5}" ]]; then
      printf '  checksum=md5=%s\n' "${expected_md5}"
      printf '  check-integrity=true\n'
    fi
  } > "${tmp_input}"

  local max_concurrent=1
  if (( ARIA2_J_CAP > 1 )); then
    max_concurrent=1
  fi

  aria2c \
    --input-file="${tmp_input}" \
    --allow-overwrite=true \
    --auto-file-renaming=false \
    --continue=true \
    --file-allocation=none \
    --max-concurrent-downloads="${max_concurrent}" \
    --split="${ARIA2_SPLIT}" \
    --max-connection-per-server="${ARIA2_MAX_CONN}" \
    --timeout="${ARIA2_TIMEOUT}" \
    --connect-timeout="${ARIA2_CONNECT_TIMEOUT}" \
    --retry-wait="${ARIA2_RETRY_WAIT}" \
    --max-tries="${ARIA2_MAX_TRIES}" \
    --summary-interval=0 \
    --console-log-level=warn

  rm -f "${tmp_input}" 2>/dev/null || true

  if [[ -n "${expected_md5}" ]]; then
    md5_ok "${out_file}" "${expected_md5}"
  else
    [[ -s "${out_file}" ]]
  fi
}

download_many_with_aria2() {
  # Downloads multiple (url, filename, md5) triples.
  # Uses aria2 input file so paired-end files can download concurrently.
  #
  # Globals expected:
  #   URLS[], OUTS[], MD5S[]
  local n="${#URLS[@]}"
  if (( n == 0 )); then
    return 1
  fi

  local tmp_input
  tmp_input="$(mktemp aria2_batch.XXXXXX.txt)"

  for (( i=0; i<n; i++ )); do
    {
      printf '%s\n' "${URLS[$i]}"
      printf '  out=%s\n' "${OUTS[$i]}"
      if [[ -n "${MD5S[$i]}" ]]; then
        printf '  checksum=md5=%s\n' "${MD5S[$i]}"
        printf '  check-integrity=true\n'
      fi
    } >> "${tmp_input}"
  done

  local max_concurrent="${n}"
  if (( max_concurrent > ARIA2_J_CAP )); then
    max_concurrent="${ARIA2_J_CAP}"
  fi

  aria2c \
    --input-file="${tmp_input}" \
    --allow-overwrite=true \
    --auto-file-renaming=false \
    --continue=true \
    --file-allocation=none \
    --max-concurrent-downloads="${max_concurrent}" \
    --split="${ARIA2_SPLIT}" \
    --max-connection-per-server="${ARIA2_MAX_CONN}" \
    --timeout="${ARIA2_TIMEOUT}" \
    --connect-timeout="${ARIA2_CONNECT_TIMEOUT}" \
    --retry-wait="${ARIA2_RETRY_WAIT}" \
    --max-tries="${ARIA2_MAX_TRIES}" \
    --summary-interval=0 \
    --console-log-level=warn

  rm -f "${tmp_input}" 2>/dev/null || true

  # Verify md5 explicitly (parallel over files if possible)
  local ok=0
  ok=0
  for (( i=0; i<n; i++ )); do
    if [[ -n "${MD5S[$i]}" ]]; then
      if ! md5_ok "${OUTS[$i]}" "${MD5S[$i]}"; then
        ok=1
      fi
    else
      if [[ ! -s "${OUTS[$i]}" ]]; then
        ok=1
      fi
    fi
  done

  return "${ok}"
}

bz2_to_gz_stream() {
  # Args: in.bz2 out.gz total_cpus
  local in_bz2="$1"
  local out_gz="$2"
  local total_cpus="$3"

  local pbz_threads=1
  local pigz_threads=1
  if (( total_cpus >= 2 )); then
    pbz_threads=$(( total_cpus / 2 ))
    pigz_threads=$(( total_cpus - pbz_threads ))
    if (( pbz_threads < 1 )); then pbz_threads=1; fi
    if (( pigz_threads < 1 )); then pigz_threads=1; fi
  fi

  pbzip2 -dc -p "${pbz_threads}" "${in_bz2}" | pigz -p "${pigz_threads}" > "${out_gz}"
}

fetch_ena_filereport() {
  # Writes ena_filereport.tsv
  local url="https://www.ebi.ac.uk/ena/portal/api/filereport"
  local fields="run_accession,fastq_ftp,fastq_md5,fastq_bytes,submitted_ftp,submitted_md5,submitted_bytes,sra_ftp,sra_md5,sra_bytes"
  curl -fsSL \
    --retry 3 \
    --retry-delay 2 \
    --connect-timeout 20 \
    --max-time 60 \
    "${url}?accession=${srr}&result=read_run&fields=${fields}&format=tsv" \
    > ena_filereport.tsv
}

parse_ena_fields() {
  # Outputs globals:
  # FASTQ_FTP FASTQ_MD5 SUB_FTP SUB_MD5 SRA_FTP SRA_MD5
  FASTQ_FTP="$(awk -F'\t' 'NR==2 {print $2}' ena_filereport.tsv || true)"
  FASTQ_MD5="$(awk -F'\t' 'NR==2 {print $3}' ena_filereport.tsv || true)"
  SUB_FTP="$(awk -F'\t' 'NR==2 {print $5}' ena_filereport.tsv || true)"
  SUB_MD5="$(awk -F'\t' 'NR==2 {print $6}' ena_filereport.tsv || true)"
  SRA_FTP="$(awk -F'\t' 'NR==2 {print $8}' ena_filereport.tsv || true)"
  SRA_MD5="$(awk -F'\t' 'NR==2 {print $9}' ena_filereport.tsv || true)"
}

try_ena_fastq_or_submitted() {
  # Builds URLS/OUTS/MD5S arrays from ENA fastq_ftp first; else submitted_ftp filtered to fastq.
  URLS=()
  OUTS=()
  MD5S=()

  local ftp_list=""
  local md5_list=""

  if [[ -n "${FASTQ_FTP}" && -n "${FASTQ_MD5}" ]]; then
    ftp_list="${FASTQ_FTP}"
    md5_list="${FASTQ_MD5}"
  elif [[ -n "${SUB_FTP}" && -n "${SUB_MD5}" ]]; then
    ftp_list="${SUB_FTP}"
    md5_list="${SUB_MD5}"
  else
    return 1
  fi

  local -a ftp_arr=()
  local -a md5_arr=()
  IFS=';' read -r -a ftp_arr <<< "${ftp_list}"
  IFS=';' read -r -a md5_arr <<< "${md5_list}"

  local n="${#ftp_arr[@]}"
  if (( n == 0 )); then
    return 1
  fi

  for (( i=0; i<n; i++ )); do
    local raw="${ftp_arr[$i]}"
    local md5=""
    if (( i < ${#md5_arr[@]} )); then
      md5="${md5_arr[$i]}"
    fi

    raw="$(printf '%s' "${raw}" | sed 's/^[[:space:]]*//; s/[[:space:]]*$//')"
    md5="$(printf '%s' "${md5}" | sed 's/^[[:space:]]*//; s/[[:space:]]*$//')"

    if [[ -z "${raw}" ]]; then
      continue
    fi

    local fname
    fname="$(basename "${raw}")"

    # If we are using submitted_ftp, filter to FASTQ-like names.
    if [[ "${ftp_list}" == "${SUB_FTP}" ]]; then
      if ! is_fastq_name "${fname}"; then
        continue
      fi
    fi

    URLS+=("$(normalise_ena_url "${raw}")")
    OUTS+=("${fname}")
    MD5S+=("${md5}")
  done

  if (( ${#URLS[@]} == 0 )); then
    return 1
  fi

  download_many_with_aria2
}

convert_bz2_fastqs_if_any() {
  shopt -s nullglob
  local -a bz2s=( ./*.fastq.bz2 ./*.fq.bz2 )
  shopt -u nullglob

  if (( ${#bz2s[@]} == 0 )); then
    return 0
  fi

  log "Converting ${#bz2s[@]} .bz2 FASTQ file(s) to .gz"
  local f
  for f in "${bz2s[@]}"; do
    local out="${f%.bz2}.gz"
    bz2_to_gz_stream "${f}" "${out}" "${cpus}"
    rm -f "${f}" 2>/dev/null || true
  done
}

try_ena_sra_then_fasterq_dump() {
  if [[ -z "${SRA_FTP}" ]]; then
    return 1
  fi

  local raw
  raw="${SRA_FTP%%;*}"
  raw="$(printf '%s' "${raw}" | sed 's/^[[:space:]]*//; s/[[:space:]]*$//')"
  if [[ -z "${raw}" ]]; then
    return 1
  fi

  local url
  url="$(normalise_ena_url "${raw}")"

  local sra_file
  sra_file="$(basename "${raw}")"
  if [[ "${sra_file}" != *.sra ]]; then
    # Still download it, but name it predictably.
    sra_file="${srr}.sra"
  fi

  log "Downloading ENA SRA file ${sra_file}"
  if ! aria2_download_with_md5 "${url}" "${sra_file}" "${SRA_MD5}"; then
    rm -f "${sra_file}" 2>/dev/null || true
    return 1
  fi

  mkdir -p tmp_srr

  log "Converting SRA -> FASTQ via fasterq-dump"
  if ! fasterq-dump \
      -e "${cpus}" \
      -t tmp_srr \
      -p \
      --outdir . \
      "${sra_file}"
  then
    rm -rf tmp_srr 2>/dev/null || true
    return 1
  fi

  rm -rf tmp_srr 2>/dev/null || true

  shopt -s nullglob
  local -a fq=( ./*.fastq ./*.fq )
  shopt -u nullglob
  if (( ${#fq[@]} > 0 )); then
    log "Compressing FASTQ with pigz"
    pigz -p "${cpus}" "${fq[@]}"
  fi

  return 0
}

try_ncbi_prefetch_then_fasterq_dump() {
  mkdir -p tmp_srr

  log "Falling back to NCBI prefetch"
  if ! prefetch \
      --output-directory . \
      --max-size u \
      "${srr}"
  then
    rm -rf tmp_srr 2>/dev/null || true
    return 1
  fi

  local sra_path=""
  if [[ -f "${srr}.sra" ]]; then
    sra_path="${srr}.sra"
  else
    sra_path="$(find . -maxdepth 3 -type f -name "${srr}.sra" -print -quit || true)"
  fi

  if [[ -z "${sra_path}" ]]; then
    rm -rf tmp_srr 2>/dev/null || true
    return 1
  fi

  # If ENA provided sra_md5, validate the downloaded .sra against it.
  if [[ -n "${SRA_MD5}" ]]; then
    log "Validating .sra against ENA sra_md5"
    if ! md5_ok "${sra_path}" "${SRA_MD5}"; then
      rm -rf tmp_srr 2>/dev/null || true
      return 1
    fi
  else
    # Otherwise do a structural validate (not md5, but better than nothing).
    vdb-validate "${sra_path}" >/dev/null 2>&1 || true
  fi

  log "Converting SRA -> FASTQ via fasterq-dump"
  if ! fasterq-dump \
      -e "${cpus}" \
      -t tmp_srr \
      -p \
      --outdir . \
      "${sra_path}"
  then
    rm -rf tmp_srr 2>/dev/null || true
    return 1
  fi

  rm -rf tmp_srr 2>/dev/null || true

  shopt -s nullglob
  local -a fq=( ./*.fastq ./*.fq )
  shopt -u nullglob
  if (( ${#fq[@]} > 0 )); then
    log "Compressing FASTQ with pigz"
    pigz -p "${cpus}" "${fq[@]}"
  fi

  return 0
}

# ---------------------- Main control flow ----------------------

# 1) Try to fetch ENA file report (small + retry).
ENA_OK=0
for (( t=1; t<=INTERNAL_RETRIES; t++ )); do
  if fetch_ena_filereport; then
    ENA_OK=1
    break
  fi
  log "ENA filereport fetch failed (try ${t}/${INTERNAL_RETRIES})"
  backoff_sleep "${t}"
done

FASTQ_FTP=""
FASTQ_MD5=""
SUB_FTP=""
SUB_MD5=""
SRA_FTP=""
SRA_MD5=""

if (( ENA_OK == 1 )); then
  parse_ena_fields
fi

# 2) Download ENA FASTQ/submitted FASTQ if available.
DONE=0
if (( ENA_OK == 1 )); then
  for (( t=1; t<=INTERNAL_RETRIES; t++ )); do
    if try_ena_fastq_or_submitted; then
      convert_bz2_fastqs_if_any
      DONE=1
      break
    fi
    log "ENA FASTQ/submitted download failed (try ${t}/${INTERNAL_RETRIES})"
    cleanup_partial
    backoff_sleep "${t}"
  done
fi

# 3) If not done, try ENA SRA -> fasterq-dump
if (( DONE == 0 )) && (( ENA_OK == 1 )); then
  for (( t=1; t<=INTERNAL_RETRIES; t++ )); do
    if try_ena_sra_then_fasterq_dump; then
      DONE=1
      break
    fi
    log "ENA SRA fallback failed (try ${t}/${INTERNAL_RETRIES})"
    cleanup_partial
    backoff_sleep "${t}"
  done
fi

# 4) If still not done, NCBI prefetch -> fasterq-dump
if (( DONE == 0 )); then
  for (( t=1; t<=INTERNAL_RETRIES; t++ )); do
    if try_ncbi_prefetch_then_fasterq_dump; then
      DONE=1
      break
    fi
    log "NCBI prefetch fallback failed (try ${t}/${INTERNAL_RETRIES})"
    cleanup_partial
    backoff_sleep "${t}"
  done
fi

if (( DONE == 0 )); then
  if (( attempt < max_retries )); then
    cleanup_partial
    exit 1
  fi
  echo "Fastq: download/conversion failed (ENA+SRA-toolkit fallback exhausted)" > FAIL.note
  cleanup_partial
  rm -f ./*.f*q* "${srr}.sra" 2>/dev/null || true
  exit 0
fi

# PacBio assembler check (only if platform is PACBIO_SMRT and assembler is missing/unknown)
final_asm="${assembler}"
if [[ "${platform}" == "PACBIO_SMRT" && ( -z "${assembler}" || "${assembler}" == "unknown" ) ]]; then
  log "Checking PacBio reads to determine assembler"
  if zcat -f ./*.f*q* 2>/dev/null \
    | awk 'NR%4==1{ h=tolower($0); if (h ~ /\/ccs([[:space:]]|$)/) { found=1; exit } } END{ exit(!found) }'
  then
    final_asm="hifi"
  else
    final_asm="pacbio"
  fi
  printf '%s\n' "${final_asm}" > assembler.txt
else
  # Always provide an assembler.txt for downstream mapping
  : > assembler.txt
fi
