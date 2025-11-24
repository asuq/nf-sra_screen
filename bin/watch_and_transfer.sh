#!/usr/bin/env bash
#
# Monitor a Nextflow run directory and, for each newly finished sample:
#   - submit an sbatch job to rsync output/$sra/$srr to DEST_DIR/$sra/$srr
#   - in that job, remove the source directory (output/$sra/$srr) after a successful transfer
#   - verify via sacct that the transfer completed successfully
#   - record the sample as processed in RUN_DIR/.processed_summary.tsv
#
# Purpose:
#   Free up space in the Nextflow run directory by moving final outputs
#   to a separate storage location and removing the per-sample output
#   directory only after a verified successful transfer.
#
# Usage:
#   watch_and_transfer.sh RUN_DIR DEST_DIR INTERVAL_MINS
#
#   RUN_DIR        Nextflow run directory (has output/summary.tsv)
#   DEST_DIR       Destination base directory for final outputs
#   INTERVAL_MINS  Positive integer; number of minutes between scans
#
# The script is intended to run inside a long-lived tmux session.

set -o errexit
set -o pipefail
set -o nounset

# ----------------------------------------------------------------------
# Logging and usage helpers
# ----------------------------------------------------------------------

log() {
    # Timestamped log messages to stderr
    printf '[%s]' "$(date +'%F %T')" >&2
    printf ' %s' "$@" >&2
    printf '\n' >&2
}

usage() {
    cat <<EOF
Usage: $(basename "$0") RUN_DIR DEST_DIR INTERVAL_MINS

  RUN_DIR        Nextflow run directory (has output/summary.tsv)
  DEST_DIR       Destination base directory
  INTERVAL_MINS  Positive integer number of minutes between checks

Internal state:
  RUN_DIR/.processed_summary.tsv : processed (sra, srr)
  RUN_DIR/.pending_copy_jobs.tsv : pending (sra, srr, job_id)

Environment:
  SACCT_CHUNK       Number of JobIDs per sacct query (default: 50)
EOF
}

# ----------------------------------------------------------------------
# Argument handling
# ----------------------------------------------------------------------

if [ "$#" -ne 3 ]; then
    usage
    exit 1
fi

RUN_DIR=$1
DEST_DIR=$2
INTERVAL_MINS=$3

# Strip trailing slashes for neat path joins
RUN_DIR=${RUN_DIR%/}
DEST_DIR=${DEST_DIR%/}

SUMMARY_FILE="$RUN_DIR/output/summary.tsv"
STATE_FILE="$RUN_DIR/.processed_summary.tsv"
PENDING_FILE="$RUN_DIR/.pending_copy_jobs.tsv"
LOCK_FILE="$RUN_DIR/.watch_and_transfer.lock"

# Validate INTERVAL_MINS as positive integer
case "$INTERVAL_MINS" in
    ''|*[!0-9]*)
        log "INTERVAL_MINS must be a positive integer (minutes), got '$INTERVAL_MINS'"
        exit 1
        ;;
esac

if [ "$INTERVAL_MINS" -le 0 ]; then
    log "INTERVAL_MINS must be greater than zero, got '$INTERVAL_MINS'"
    exit 1
fi

SLEEP_SECONDS=$((INTERVAL_MINS * 60))

# Basic sanity checks
if [ ! -d "$RUN_DIR" ]; then
    log "RUN_DIR does not exist: $RUN_DIR"
    exit 1
fi

if [ ! -d "$DEST_DIR" ]; then
    log "DEST_DIR does not exist: $DEST_DIR"
    exit 1
fi

for cmd in sbatch sacct flock rsync; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        log "$cmd not found in PATH; this script requires $cmd."
        exit 1
    fi
done

# Acquire an exclusive, non-blocking lock on the run directory.
# The lock is held on file descriptor 9 for the lifetime of this process.
exec 9>"$LOCK_FILE" || {
    log "Unable to open lock file $LOCK_FILE"
    exit 1
}

if ! flock -n 9; then
    log "Another instance of this watcher appears to be running for $RUN_DIR; exiting."
    exit 1
fi

# Ensure state files exist
if [ ! -e "$STATE_FILE" ]; then
    : > "$STATE_FILE"
    log "Initialised state file: $STATE_FILE"
fi

if [ ! -e "$PENDING_FILE" ]; then
    : > "$PENDING_FILE"
    log "Initialised pending jobs file: $PENDING_FILE"
fi

# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------

is_in_pending() {
    # Return 0 if (sra, srr) exists in PENDING_FILE, using strict tab-field match.
    local sra="$1" srr="$2"
    [ -s "$PENDING_FILE" ] || return 1
    awk -F'\t' -v sra="$sra" -v srr="$srr" '
        ($1 == sra && $2 == srr) { found=1; exit }
        END { exit (found ? 0 : 1) }
    ' "$PENDING_FILE"
}

check_pending_jobs() {
    # Reconcile all pending (sra, srr, job_id) entries with batched sacct queries.
    # Successful jobs append (sra, srr) to STATE_FILE and are removed from pending.
    if [ ! -r "$PENDING_FILE" ] || [ ! -s "$PENDING_FILE" ]; then
        return
    fi

    # Collect job IDs into an array
    local -a job_ids=()
    while IFS=$'\t' read -r sra srr job_id; do
        [ -z "${sra:-}" ] && continue
        [ -z "${srr:-}" ] && continue
        [ -z "${job_id:-}" ] && continue
        job_ids+=("$job_id")
    done < "$PENDING_FILE"

    # Batched sacct query in chunks (default 50 IDs per call; override with SACCT_CHUNK)
    local sacct_output=""
    local chunk_size_raw="${SACCT_CHUNK:-50}"
		case "$chunk_size_raw" in
				''|*[!0-9]*)
						log "SACCT_CHUNK must be a positive integer; got '$chunk_size_raw'. Falling back to 50."
						chunk_size_raw=50
						;;
		esac
		if [ "$chunk_size_raw" -le 0 ]; then
				log "SACCT_CHUNK must be > 0; got '$chunk_size_raw'. Falling back to 50."
				chunk_size_raw=50
		fi
		local chunk_size="$chunk_size_raw"


    local total="${#job_ids[@]}"
    local i=0
    while [ "$i" -lt "$total" ]; do
        local -a chunk=( "${job_ids[@]:$i:$chunk_size}" )
        # Join chunk into comma-separated list
        local csv
        csv="$(printf '%s,' "${chunk[@]}")"; csv="${csv%,}"
        # Query sacct; append any output. Use '|| true' to avoid exiting on transient failures.
        local out
        out="$(sacct -j "$csv" --format=JobIDRaw,State,ExitCode --noheader --parsable2 2>/dev/null || true)"
        if [ -n "$out" ]; then
            if [ -n "$sacct_output" ]; then
                sacct_output="${sacct_output}"$'\n'"${out}"
            else
                sacct_output="${out}"
            fi
        fi
        i=$(( i + chunk_size ))
    done

    # Build mapping: base_jobid -> "base_jobid|State|ExitCode" (prefer exact base row)
    local map_file
    map_file="$(mktemp "${PENDING_FILE}.sacctmap.XXXXXX")"
    if [ -n "$sacct_output" ]; then
        printf '%s\n' "$sacct_output" | awk -F'|' '
            {
                raw=$1; st=$2; ex=$3;
                base=raw; sub(/\..*$/, "", base);
                rank=(raw==base)?2:1
                if (!(base in best) || rank>best[base]) {
                    best[base]=rank; state[base]=st; exitc[base]=ex
                }
            }key
            END {
                for (b in best) print b "|" state[b] "|" exitc[b]
            }
        ' > "$map_file"
    else
        : > "$map_file"
    fi

    # New pending file (atomic mv)
    local tmp_pending
    tmp_pending="$(mktemp "${PENDING_FILE}.tmp.XXXXXX")"
    : > "$tmp_pending"

    # Second pass: reconcile rows
    while IFS=$'\t' read -r sra srr job_id; do
        [ -z "${sra:-}" ] && continue
        [ -z "${srr:-}" ] && continue
        [ -z "${job_id:-}" ] && continue

        line="$(awk -F'|' -v id="$job_id" '$1==id {print; exit}' "$map_file" || true)"

        if [ -z "$line" ]; then
            # No accounting record yet; keep pending
            printf '%s\t%s\t%s\n' "$sra" "$srr" "$job_id" >> "$tmp_pending"
            continue
        fi

        job_state="$(printf '%s\n' "$line" | awk -F'|' '{print $2}')"
        job_exit="$(printf '%s\n' "$line" | awk -F'|' '{print $3}')"

        case "$job_state" in
            PENDING|CONFIGURING|REQUEUED|RESIZING|RUNNING|COMPLETING|SUSPENDED)
                # Still in progress; keep pending
                printf '%s\t%s\t%s\n' "$sra" "$srr" "$job_id" >> "$tmp_pending"
                ;;
            COMPLETED)
                case "$job_exit" in
                    0:*)
                        log "Transfer job $job_id for $sra/$srr completed successfully (ExitCode=$job_exit)."

                        # Clean up the Slurm log for a successful transfer
                        log_file="${RUN_DIR}/slurm-${job_id}.out"
                        if [ -f "$log_file" ]; then
                            rm -vf -- "$log_file" || \
                              log "Warning: failed to remove log file $log_file"
                        fi

                        # Record the sample as fully processed so we never handle it again.
                        printf '%s\t%s\n' "$sra" "$srr" >> "$STATE_FILE"
                        ;;
                    *)
                        log "Transfer job $job_id for $sra/$srr completed with non-zero ExitCode=$job_exit; will be resubmitted in a later cycle."
                        ;;
                esac
                ;;
            *)
                # Any other terminal state (FAILED, CANCELLED, TIMEOUT, etc.)
                log "Transfer job $job_id for $sra/$srr reached terminal State=$job_state ExitCode=$job_exit; will be resubmitted in a later cycle."
                ;;
        esac
    done < "$PENDING_FILE"

    mv "$tmp_pending" "$PENDING_FILE"
    rm -f "$map_file"
}

# ----------------------------------------------------------------------
# Core worker: process any new samples found in summary.tsv
# ----------------------------------------------------------------------
move_output_to_storage() {
    # For each (sra, srr) not yet in STATE_FILE and not already pending:
    #   - if note starts with "did not match the criteria": mark processed (no transfer)
    #   - else if note is non-empty and no output dir exists: mark processed (no transfer)
    #   - else, ensure output directory exists, then submit an sbatch rsync job
    #     that transfers the directory and removes the source on success.

    if [ ! -r "$SUMMARY_FILE" ]; then
        log "summary.tsv not readable: $SUMMARY_FILE -- skipping this cycle."
        return
    fi

    awk -F'\t' 'NR>1 { OFS="\t"; print $1,$2,$8 }' "$SUMMARY_FILE" | \
    while IFS=$'\t' read -r sra srr note; do
        # Skip empty or malformed lines
        if [ -z "${sra:-}" ] || [ -z "${srr:-}" ]; then
            continue
        fi

        key="${sra}"$'\t'"${srr}"

        # Already processed?
        if grep -Fqx "$key" "$STATE_FILE"; then
            continue
        fi

        # Already pending?
        if is_in_pending "$sra" "$srr"; then
            continue
        fi

        # Normalise note (trim leading/trailing whitespace)
        note_trimmed="$(printf '%s' "${note:-}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"

        sample_src="$RUN_DIR/output/$sra/$srr"

        # 1. Special case: "did not match the criteria: ..."
        # These are from the initial metadata filter. They are logically
        # finished, but there is no downstream analysis to transfer.
        case "$note_trimmed" in
            "did not match the criteria"*)
                log "Sample $sra/$srr did not pass the selection criteria:"
                log "  note      : $note_trimmed"
                log "  sample_src: $sample_src"
                log "Treating as filtered sample; marking as processed (no transfer job submitted)."

                printf '%s\t%s\n' "$sra" "$srr" >> "$STATE_FILE"
                continue
                ;;
        esac

        # 2. Any other note + NO output directory
        # e.g. LOG_FAILED_PROCESS with no outputs. Log it and mark processed.
        if [ -n "$note_trimmed" ] && [ ! -d "$sample_src" ]; then
            log "Sample $sra/$srr has summary note but no output directory:"
            log "  note      : $note_trimmed"
            log "  sample_src: $sample_src"
            log "Treating as filtered/failed with no outputs; marking as processed."

            printf '%s\t%s\n' "$sra" "$srr" >> "$STATE_FILE"
            continue
        fi

        # 3. Normal completed sample (with or without note)
        log "Handling finished sample: sra=$sra  srr=$srr"

        # If we still have no output directory and no note, this is a weird partial state.
        if [ ! -d "$sample_src" ]; then
					if [ -d "$sample_dest" ]; then
						log "Output dir missing but destination exists for $sra/$srr:"
						log "  sample_src : $sample_src"
						log "  sample_dest: $sample_dest"
						log "Assuming prior successful transfer; marking as processed."
						printf '%s\t%s\n' "$sra" "$srr" >> "$STATE_FILE"
					else
            log "Output directory missing for $sra/$srr: $sample_src"
            log "Sample will not be marked as done; will retry next cycle."
					fi
          continue
        fi

        # If there is a note *and* an output directory, we keep the outputs
        # but mention the note for traceability.
        if [ -n "$note_trimmed" ]; then
            log "Sample $sra/$srr has note in summary but also an output directory:"
            log "  note=$note_trimmed"
        fi

        sample_dest_parent="$DEST_DIR/$sra"
        sample_dest="$sample_dest_parent/$srr"

        # Transfer: rsync -a (copy) then rm -rf (remove source) if rsync succeeds.
        transfer_cmd="mkdir -p \"$sample_dest\" && rsync -a \"$sample_src\"/ \"$sample_dest\"/ && rm -rf -- \"$sample_src\""

        log "Submitting transfer job for $sra/$srr:"
        log "  from: $sample_src"
        log "  to  : $sample_dest"

        # Capture output safely under errexit (do not exit on sbatch failure)
        sbatch_output="$(sbatch -p datacp --parsable --job-name="move_${srr}" --wrap="$transfer_cmd" 2>&1 || true)"
        if [ -z "$sbatch_output" ]; then
            log "sbatch submission produced no output for $sra/$srr"
            log "Sample will not be marked as done; will retry next cycle."
            continue
        fi

        job_id="$(printf '%s\n' "$sbatch_output" | awk -F';' 'NR==1 {print $1}' | tr -d '[:space:]')"
        if [ -z "$job_id" ] || ! printf '%s' "$job_id" | grep -Eq '^[0-9]+'; then
            log "Could not parse job ID from sbatch output for $sra/$srr: $sbatch_output"
            log "Sample will not be marked as done; will retry next cycle."
            continue
        fi

        log "Submitted transfer job $job_id for $sra/$srr; it will be monitored in subsequent cycles."
        printf '%s\t%s\t%s\n' "$sra" "$srr" "$job_id" >> "$PENDING_FILE"
    done
}

# ----------------------------------------------------------------------
# Main loop
# ----------------------------------------------------------------------
log "Starting Nextflow watcher."
log "  Run dir : $RUN_DIR"
log "  Dest dir: $DEST_DIR"
log "  Interval: $INTERVAL_MINS minute(s) (= $SLEEP_SECONDS seconds)"
log "  State   : $STATE_FILE"
log "  Pending : $PENDING_FILE"

while :; do
    check_pending_jobs
    move_output_to_storage
    log "Cycle complete; sleeping for $INTERVAL_MINS minute(s)."
    sleep "$SLEEP_SECONDS"
done
