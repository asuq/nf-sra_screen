#!/usr/bin/env bash
# Usage:
#   run_singlem.sh \
#       --sandpiper-decision PASS|RUN_SINGLEM \
#       --valid-taxa validated_taxa.csv \
#       --singlem-db /path/to/singlem/metapackage \
#       --cpus N \
#       --attempt A \
#       --max-retries M \
#       --reads READS.fastq.gz [...]
#
# Produces:
#   - singlem_taxonomic_profile.tsv
#   - singlem_taxonomic_profile_krona*
#   - singlem_output.tsv
#   - reads_ok/*.f*q*
#   - FAIL.note on terminal failure or logical negative

set -euo pipefail

sandpiper_dec=""
valid_taxa=""
singlem_db=""
cpus=1
attempt=0
max_retries=1
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
    --sandpiper-decision) sandpiper_dec="$2"; shift 2 ;;
    --valid-taxa)         valid_taxa="$2"; shift 2 ;;
    --singlem-db)         singlem_db="$2"; shift 2 ;;
    --cpus)               cpus="$2"; shift 2 ;;
    --attempt)            attempt="$2"; shift 2 ;;
    --max-retries)        max_retries="$2"; shift 2 ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1 ;;
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

# Skip SingleM entirely if Sandpiper already passed
if [[ "$sandpiper_dec" == "PASS" ]]; then
  echo "Reads already passed Sandpiper check; skipping SingleM"
  mkdir -p reads_ok
  for f in "${read_files[@]}"; do
    # Keep behaviour: real copies, not symlinks
    cp -v "$f" "reads_ok/$(basename "$f")"
  done
  exit 0
fi

if [[ ${#read_files[@]} -eq 0 || -z "$valid_taxa" || -z "$singlem_db" ]]; then
  fail "SingleM: missing required arguments (--reads, --valid-taxa, --singlem-db)"
fi

# Build phyla_to_check.txt from validated_taxa.csv
make_phyla_list_from_validated_taxa.sh "$valid_taxa" phyla_to_check.txt \
  || fail "SingleM: failed to parse validated_taxa"

# Parse reads: paired‑end if ≥2, else single‑end
if [[ "${#read_files[@]}" -ge 2 ]]; then
  R1="${read_files[0]}"
  R2="${read_files[1]}"
else
  R1="${read_files[0]}"
  R2=""
fi

# Run SingleM pipe
if [[ -n "$R1" && -n "$R2" ]]; then
  singlem pipe \
    -1 "$R1" \
    -2 "$R2" \
    --taxonomic-profile singlem_taxonomic_profile.tsv \
    --taxonomic-profile-krona singlem_taxonomic_profile_krona \
    --metapackage "$singlem_db" \
    --threads "$cpus" \
    || fail "SingleM: pipe failed"
else
  singlem pipe \
    -1 "$R1" \
    --taxonomic-profile singlem_taxonomic_profile.tsv \
    --taxonomic-profile-krona singlem_taxonomic_profile_krona \
    --metapackage "$singlem_db" \
    --threads "$cpus" \
    || fail "SingleM: pipe failed"
fi

if [[ ! -s singlem_taxonomic_profile.tsv ]]; then
  fail "SingleM: empty taxonomic profile"
fi

# SingleM summarise
singlem summarise \
  --input-taxonomic-profile singlem_taxonomic_profile.tsv \
  --output-species-by-site-relative-abundance-prefix singlem_taxonomic_profile_summary \
  || fail "SingleM: summarise failed"

# ---------------------- Inline phylum check logic ---------------------- #

rc=0
if check_singlem_phyla.py \
      -i "singlem_taxonomic_profile.tsv" \
      -p "phyla_to_check.txt" \
      -o "singlem_output.tsv"; then
  rc=0
else
  rc=$?
fi

case "$rc" in
  0)
		shopt -s nullglob
		mkdir -p reads_ok
		for f in *.f*q*; do
			ln -sf "../$f" "reads_ok/$(basename "$f")"
		done
		shopt -u nullglob
    ;;
  1)
    fail "SingleM: phylum check internal error"
    ;;
  2)
    echo "SingleM: No target phyla detected" > FAIL.note
		exit 0
    ;;
	*)
    fail "SingleM: phylum check unexpected exit code ${rc}"
    ;;
esac