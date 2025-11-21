#!/usr/bin/env bash
#
# Monitor a Nextflow run directory and, for each newly finished sample:
#   - delete all work/ directories related to that sample (based on .nextflow.log)
#   - transfer output/$sra/$srr to DEST_DIR using rsync -a; rm -rf via sbatch on datacp partition
#   - verify via sacct that the transfer completed successfully
#   - record the sample as processed in RUN_DIR/.processed_summary.tsv
#
# Purpose:
#   Free up space in the Nextflow run directory by deleting work/ data
#   for samples that have completed, and move their final outputs to
#   a separate storage location.
#
# Usage:
#   watch_and_move.sh RUN_DIR DEST_DIR INTERVAL_MINS
#
#   RUN_DIR        Nextflow run directory (has .nextflow.log, work/, output/summary.tsv)
#   DEST_DIR       Destination base directory for final outputs
#   INTERVAL_MINS Positive integer; number of minutes between scans
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

  RUN_DIR        Nextflow run directory (has .nextflow.log, work/, output/summary.tsv)
  DEST_DIR       Destination base directory
  INTERVAL_MINS  Positive integer number of minutes between checks

Internal state:
  RUN_DIR/.processed_summary.tsv : processed (sra, srr)
  RUN_DIR/.pending_copy_jobs.tsv : pending (sra, srr, job_id)
  RUN_DIR/.workdirs_index.tsv    : per-cycle index mapping (sra, srr) -> work_dir

Environment:
  SACCT_CHUNK   Number of JobIDs per sacct query (default: 50)
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
NEXTFLOW_LOG="$RUN_DIR/.nextflow.log"
WORK_DIR="$RUN_DIR/work"
STATE_FILE="$RUN_DIR/.processed_summary.tsv"
PENDING_FILE="$RUN_DIR/.pending_copy_jobs.tsv"
LOCK_FILE="$RUN_DIR/.watch_and_move.lock"
WORK_INDEX_FILE="$RUN_DIR/.workdirs_index.tsv"

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

if [ ! -d "$WORK_DIR" ]; then
    log "work/ directory does not exist under RUN_DIR: $WORK_DIR"
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
build_workdir_index() {
	# Parse .nextflow.log once per cycle and build an index mapping (sra, srr) -> work_dir.
	# Output: RUN_DIR/.workdirs_index.tsv (tab-separated: sra, srr, work_dir).
    if [ ! -r "$NEXTFLOW_LOG" ]; then
        : > "$WORK_INDEX_FILE"
        log ".nextflow.log not readable; workdir index will be empty this cycle."
        return
    fi

    # Create temp file in the same directory for atomic mv
    local tmp_index
    tmp_index="$(mktemp "${WORK_INDEX_FILE}.tmp.XXXXXX")"

    # Extract lines with "(SRA:SRR)" and a 'workDir:' token; print sra \t srr \t work_path
    LC_ALL=C awk '
        {
            if (match($0, /\(([[:alnum:]_.:-]+):([[:alnum:]_.:-]+)\)/, m)) {
                sra = m[1]; srr = m[2];
                for (i = 1; i <= NF; i++) {
                    if ($i == "workDir:") {
                        if (i + 1 <= NF) {
                            val = $(i + 1)
                            sub(/[,\]]+$/, "", val)
                            print sra "\t" srr "\t" val
                        }
                    }
                }
            }
        }
    ' "$NEXTFLOW_LOG" | sort -u > "$tmp_index"

    mv "$tmp_index" "$WORK_INDEX_FILE"
    log "Rebuilt workdir index: $WORK_INDEX_FILE"
}

is_in_pending() {
	# Return 0 if (sra, srr) exists in PENDING_FILE, using strict tab-field match.
    local sra="$1" srr="$2"
    [ -s "$PENDING_FILE" ] || return 1
    awk -F'\t' -v sra="$sra" -v srr="$srr" '
        ($1 == sra && $2 == srr) { found=1; exit }
        END { exit (found ? 0 : 1) }
    ' "$PENDING_FILE"
}

delete_sample_work_dirs() {
	#  Use the per-cycle workdir index to find work dirs for (sra, srr) and delete only those under $WORK_DIR/.
    local sra="$1" srr="$2"

    if [ ! -r "$WORK_INDEX_FILE" ] || [ ! -s "$WORK_INDEX_FILE" ]; then
        log "Workdir index missing/empty; skipping work/ deletion for $sra/$srr."
        return
    fi

    log "Collecting work directories for sample sra=$sra srr=$srr from index"

    awk -F'\t' -v sra="$sra" -v srr="$srr" '
        $1 == sra && $2 == srr { print $3 }
    ' "$WORK_INDEX_FILE" | sort -u | \
    while IFS= read -r work_path; do
        [ -z "$work_path" ] && continue
        case "$work_path" in
            "$WORK_DIR"/*)
                if [ -d "$work_path" ]; then
                    log "Deleting work dir for $sra/$srr: $work_path"
                    rm -rf -- "$work_path"
                fi
                ;;
            *)
								# Safety fence: never delete anything outside the run's work/
                log "Skipping suspicious work path for $sra/$srr (outside $WORK_DIR): $work_path"
                ;;
        esac
    done
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
    local chunk_size="${SACCT_CHUNK:-50}"
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
            }
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
            PENDING|CONFIGURING|RUNNING|COMPLETING|SUSPENDED)
								# Still in progress; keep pending
                printf '%s\t%s\t%s\n' "$sra" "$srr" "$job_id" >> "$tmp_pending"
                ;;
            COMPLETED)
                case "$job_exit" in
                    0:*)
                        log "Transfer job $job_id for $sra/$srr completed successfully (ExitCode=$job_exit)."
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
    #   - ensure output directory exists (or handle note-only rows)
    #   - delete all work/ dirs for that sample (using the per-cycle index)
    #   - submit an sbatch rsync job
    #   - record (sra, srr, job_id) in PENDING_FILE

    if [ ! -r "$SUMMARY_FILE" ]; then
        log "summary.tsv not readable: $SUMMARY_FILE -- skipping this cycle."
        return
    fi

    awk 'NR>1{print}' "$SUMMARY_FILE" | \
    while IFS=$'\t' read -r sra srr platform model strategy assembler counts note; do
        # Skip empty or malformed lines
        if [ -z "${sra:-}" ] || [ -z "${srr:-}" ]; then
            continue
        fi

        key=${sra}$'\t'${srr}

        # Already processed?
        if grep -Fqx "$key" "$STATE_FILE"; then
            continue
        fi

        # Already pending?
        if is_in_pending "$sra" "$srr"; then
            continue
        fi

        log "Handling finished sample: sra=$sra  srr=$srr"

        sample_src="$RUN_DIR/output/$sra/$srr"
        if [ ! -d "$sample_src" ]; then
            # If there is a note in summary.tsv but no per-SRR output directory,
            # the pipeline has already decided this SRR is "done with no output".
            # Typical case: filter_sra.sh skipped the SRR and wrote
            #   note = "did not match the criteria: <reason>"
            if [ -n "${note:-}" ]; then
                log "Output directory missing for $sra/$srr: $sample_src, but summary note='$note'."
                log "Sample will be marked as DONE with no output to transfer."
                printf '%s\t%s\n' "$sra" "$srr" >> "$STATE_FILE"
                continue
            fi

            # If note is empty and directory is missing, the pipeline hasn't
            # finished this SRR yet; keep the old behaviour and retry later.
            log "Output directory missing for $sra/$srr: $sample_src"
            log "Sample will not be marked as done; will retry next cycle."
            continue
        fi

        # Delete all work directories associated with this sample.
        delete_sample_work_dirs "$sra" "$srr"

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
log "  WorkIdx : $WORK_INDEX_FILE"

while :; do
    build_workdir_index
    check_pending_jobs
    move_output_to_storage
    log "Cycle complete; sleeping for $INTERVAL_MINS minute(s)."
    sleep "$SLEEP_SECONDS"
done
