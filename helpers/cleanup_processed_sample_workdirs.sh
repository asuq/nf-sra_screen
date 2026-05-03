#!/usr/bin/env bash
#
# Delete Nextflow work directories for samples whose trace rows all finished.
#
# The helper reads only trace.tsv. It groups SRA:SRR-style tags by sample,
# checks the final observed task status for each sample, writes the processed
# sample list, and then deletes matching work directories.

set -euo pipefail

DRY_RUN=0
FILTER_SRA=""
FILTER_SRR=""
PROCESSED_OUT=""
TRACE_FILE=""
WORK_ROOT=""
WORK_ROOT_CANONICAL=""
WORK_ROOT_EXPLICIT=0

usage() {
    # Print command-line usage.
    cat <<EOF
Usage: $(basename "$0") TRACE_TSV [options]

Options:
  --dry-run               Write processed samples and report workdirs without deleting
  --work-root WORK_DIR    Fence deletion to this work directory root
  --sra SRA               Restrict cleanup to one SRA accession
  --srr SRR               Restrict cleanup to one SRR accession
  --processed-out FILE    Processed sample TSV to write (default: beside TRACE_TSV)
  -h, --help              Show this help message

The trace file must be tab-separated and include tag, status, and workdir
headers. Rows with tags shaped SRA:SRR or SRA:SRR:* are grouped by sample.
Without --work-root, deletion is fenced under the inferred launch work/ root.
EOF
}

fail() {
    # Print an error and exit.
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

log() {
    # Print a status line to stderr.
    printf '%s\n' "$*" >&2
}

require_command() {
    # Require an executable command to be available.
    command -v "$1" >/dev/null 2>&1 || fail "required command not found: $1"
}

canonical_existing_path() {
    # Resolve an existing path to an absolute physical path.
    local path=$1
    local dir
    local base

    [ -e "$path" ] || return 1

    dir=$(dirname "$path")
    base=$(basename "$path")
    (
        cd "$dir" >/dev/null 2>&1
        printf '%s/%s\n' "$(pwd -P)" "$base"
    )
}

canonical_output_path() {
    # Resolve an output path whose parent directory already exists.
    local path=$1
    local dir
    local base

    dir=$(dirname "$path")
    base=$(basename "$path")
    [ -d "$dir" ] || fail "output directory does not exist: $dir"

    (
        cd "$dir" >/dev/null 2>&1
        printf '%s/%s\n' "$(pwd -P)" "$base"
    )
}

infer_launch_dir() {
    # Infer the launch directory used for relative trace workdir values.
    local trace_dir=$1
    local base

    base=$(basename "$trace_dir")
    if [ "$base" = "execution-reports" ]; then
        dirname "$trace_dir"
    else
        printf '%s\n' "$trace_dir"
    fi
}

normalise_trace_workdir() {
    # Print an absolute candidate workdir path from a trace workdir field.
    local workdir=$1

    case "$workdir" in
        /*)
            printf '%s\n' "$workdir"
            ;;
        *)
            printf '%s/%s\n' "$LAUNCH_DIR" "$workdir"
            ;;
    esac
}

is_safe_workdir_child() {
    # Return success when a resolved workdir is safe to delete.
    local resolved=$1
    case "$resolved" in
        "$WORK_ROOT_CANONICAL"/*)
            return 0
            ;;
        *)
            return 1
            ;;
    esac
}

parse_trace() {
    # Parse trace rows into latest task statuses and all observed workdirs.
    local final_rows=$1
    local all_workdirs=$2

    : > "$all_workdirs"

    awk -F'\t' \
        -v final_rows="$final_rows" \
        -v all_workdirs="$all_workdirs" \
        -v filter_sra="$FILTER_SRA" \
        -v filter_srr="$FILTER_SRR" '
        function upper(value) {
            return toupper(value)
        }
        function numeric(value) {
            return (value ~ /^[0-9]+$/) ? value + 0 : 0
        }
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                if ($i == "tag") tag_col = i
                if ($i == "status") status_col = i
                if ($i == "workdir") workdir_col = i
                if ($i == "attempt") attempt_col = i
                if ($i == "task_id") task_id_col = i
                if ($i == "name") name_col = i
                if ($i == "hash") hash_col = i
                if ($i == "process") process_col = i
            }
            if (!tag_col || !status_col || !workdir_col) {
                print "ERROR: trace TSV must include tag, status, and workdir headers" > "/dev/stderr"
                exit 2
            }
            next
        }
        {
            tag = $tag_col
            if (tag == "") next

            part_count = split(tag, parts, ":")
            if (part_count < 2) next

            sra = parts[1]
            srr = parts[2]
            if (sra == "" || srr == "") next
            if (filter_sra != "" && sra != filter_sra) next
            if (filter_srr != "" && srr != filter_srr) next

            sample = sra "\t" srr
            samples_seen++

            workdir = $workdir_col
            if (workdir != "") {
                print sample "\t" workdir >> all_workdirs
            }

            status = upper($status_col)
            if (name_col && $name_col != "") {
                task_key = $name_col
            }
            else if (hash_col && $hash_col != "") {
                task_key = (process_col ? $process_col : "") "|" tag "|" $hash_col
            }
            else if (attempt_col || task_id_col) {
                task_key = (process_col ? $process_col : "") "|" tag "|" workdir
            }
            else {
                task_key = sample "|" NR
            }

            score = numeric(attempt_col ? $attempt_col : "") * 1000000000 + \
                    numeric(task_id_col ? $task_id_col : "")
            if (score == 0) {
                score = NR
            }

            key = sample "\t" task_key
            if (!(key in seen)) {
                keys[++key_count] = key
                key_sample[key] = sample
                seen[key] = 1
            }
            if (!(key in best_score) || score >= best_score[key]) {
                best_score[key] = score
                best_status[key] = status
            }
        }
        END {
            if (NR == 0) {
                print "ERROR: trace TSV is empty" > "/dev/stderr"
                exit 2
            }
            for (i = 1; i <= key_count; i++) {
                key = keys[i]
                print key_sample[key] "\t" best_status[key] > final_rows
            }
        }
    ' "$TRACE_FILE"
}

write_eligible_samples() {
    # Write samples whose final observed task rows are all complete.
    local final_rows=$1
    local eligible_samples=$2

    awk -F'\t' '
        {
            sample = $1 "\t" $2
            status = $3
            if (!(sample in seen)) {
                seen[sample] = 1
                samples[++sample_count] = sample
                sra[sample] = $1
                srr[sample] = $2
            }
            if (status != "COMPLETED" && status != "CACHED") {
                blocked[sample] = 1
                if (bad_status[sample] == "") {
                    bad_status[sample] = status
                }
                else if (index("," bad_status[sample] ",", "," status ",") == 0) {
                    bad_status[sample] = bad_status[sample] "," status
                }
            }
        }
        END {
            for (i = 1; i <= sample_count; i++) {
                sample = samples[i]
                if (blocked[sample]) {
                    print "Skipping sample " sample ": final status includes " bad_status[sample] > "/dev/stderr"
                    continue
                }
                print sra[sample] "\t" srr[sample]
            }
        }
    ' "$final_rows" > "$eligible_samples"
}

filter_eligible_workdirs() {
    # Keep all observed workdirs belonging to eligible samples.
    local eligible_samples=$1
    local all_workdirs=$2
    local raw_workdirs=$3

    awk -F'\t' -v eligible_samples="$eligible_samples" '
        BEGIN {
            while ((getline line < eligible_samples) > 0) {
                split(line, fields, FS)
                if (fields[1] != "" && fields[2] != "") {
                    eligible[fields[1] FS fields[2]] = 1
                }
            }
            close(eligible_samples)
        }
        {
            sample = $1 FS $2
            if (sample in eligible && $3 != "") {
                print
            }
        }
    ' "$all_workdirs" | sort -u > "$raw_workdirs"
}

collect_safe_workdirs() {
    # Resolve, fence, and deduplicate matching work directories.
    local raw_workdirs=$1
    local safe_workdirs=$2
    local sra
    local srr
    local path
    local candidate
    local resolved

    : > "$safe_workdirs"

    while IFS=$'\t' read -r sra srr path; do
        [ -n "${sra:-}" ] || continue
        [ -n "${srr:-}" ] || continue
        [ -n "${path:-}" ] || continue

        candidate=$(normalise_trace_workdir "$path")
        if ! resolved=$(canonical_existing_path "$candidate"); then
            log "Skipping missing workdir from trace: $candidate"
            continue
        fi

        if is_safe_workdir_child "$resolved"; then
            printf '%s\t%s\t%s\n' "$sra" "$srr" "$resolved" >> "$safe_workdirs"
        else
            if [ "$WORK_ROOT_EXPLICIT" -eq 1 ]; then
                log "Skipping workdir outside --work-root: $resolved"
            else
                log "Skipping workdir outside inferred work root; pass --work-root if intended: $resolved"
            fi
        fi
    done < "$raw_workdirs"

    sort -u "$safe_workdirs" -o "$safe_workdirs"
}

write_processed_samples() {
    # Write the processed sample list before deleting any workdirs.
    local eligible_samples=$1
    local safe_workdirs=$2
    local tmp_output
    local deleted_default
    local dry_run

    deleted_default=true
    dry_run=false
    if [ "$DRY_RUN" -eq 1 ]; then
        deleted_default=false
        dry_run=true
    fi

    tmp_output=$(mktemp "${PROCESSED_OUT}.tmp.XXXXXX")
    awk -F'\t' \
        -v safe_workdirs="$safe_workdirs" \
        -v deleted_default="$deleted_default" \
        -v dry_run="$dry_run" \
        -v trace_file="$TRACE_FILE" '
        BEGIN {
            while ((getline line < safe_workdirs) > 0) {
                split(line, fields, FS)
                sample = fields[1] FS fields[2]
                if (fields[3] != "" && !(sample FS fields[3] in counted)) {
                    counted[sample FS fields[3]] = 1
                    workdir_count[sample]++
                }
            }
            close(safe_workdirs)
            print "sra\tsrr\tworkdir_count\tdeleted_default\tdry_run\ttrace_file"
        }
        {
            sample = $1 FS $2
            print $1 "\t" $2 "\t" (workdir_count[sample] + 0) "\t" deleted_default "\t" dry_run "\t" trace_file
        }
    ' "$eligible_samples" > "$tmp_output"
    mv "$tmp_output" "$PROCESSED_OUT"
}

directory_size_kib() {
    # Print a directory size estimate in KiB.
    du -sk "$1" 2>/dev/null | awk '{ print $1 }'
}

cleanup_workdirs() {
    # Delete or report the fenced work directories.
    local safe_workdirs=$1
    local unique_paths=$2
    local count=0
    local total_kib=0
    local size_kib
    local path

    cut -f 3 "$safe_workdirs" | sort -u > "$unique_paths"

    while IFS= read -r path; do
        [ -n "$path" ] || continue

        size_kib=$(directory_size_kib "$path")
        size_kib=${size_kib:-0}
        total_kib=$((total_kib + size_kib))
        count=$((count + 1))

        if [ "$DRY_RUN" -eq 1 ]; then
            printf 'DRY-RUN would delete\t%s\t%s KiB\n' "$path" "$size_kib"
        else
            printf 'Deleting\t%s\t%s KiB\n' "$path" "$size_kib"
            rm -rf -- "$path"
        fi
    done < "$unique_paths"

    if [ "$DRY_RUN" -eq 1 ]; then
        log "DRY-RUN summary: processed_samples=$(wc -l < "$PROCESSED_OUT" | awk '{print $1 - 1}') matched_workdirs=$count total_size=${total_kib} KiB"
    else
        log "Cleanup summary: processed_samples=$(wc -l < "$PROCESSED_OUT" | awk '{print $1 - 1}') deleted_workdirs=$count total_size=${total_kib} KiB"
    fi
}

while [ "$#" -gt 0 ]; do
    case "$1" in
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        --work-root)
            [ "$#" -ge 2 ] || fail "--work-root requires a value"
            WORK_ROOT=$2
            shift 2
            ;;
        --work-root=*)
            WORK_ROOT=${1#*=}
            shift
            ;;
        --sra)
            [ "$#" -ge 2 ] || fail "--sra requires a value"
            FILTER_SRA=$2
            shift 2
            ;;
        --sra=*)
            FILTER_SRA=${1#*=}
            shift
            ;;
        --srr)
            [ "$#" -ge 2 ] || fail "--srr requires a value"
            FILTER_SRR=$2
            shift 2
            ;;
        --srr=*)
            FILTER_SRR=${1#*=}
            shift
            ;;
        --processed-out)
            [ "$#" -ge 2 ] || fail "--processed-out requires a value"
            PROCESSED_OUT=$2
            shift 2
            ;;
        --processed-out=*)
            PROCESSED_OUT=${1#*=}
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        -*)
            fail "unknown option: $1"
            ;;
        *)
            if [ -n "$TRACE_FILE" ]; then
                fail "unexpected positional argument: $1"
            fi
            TRACE_FILE=$1
            shift
            ;;
    esac
done

[ -n "$TRACE_FILE" ] || { usage >&2; exit 1; }
[ -r "$TRACE_FILE" ] || fail "trace file is not readable: $TRACE_FILE"

require_command awk
require_command basename
require_command cut
require_command dirname
require_command du
require_command mktemp
require_command mv
require_command rm
require_command sort
require_command wc

TRACE_FILE=$(canonical_existing_path "$TRACE_FILE")
TRACE_DIR=$(dirname "$TRACE_FILE")
LAUNCH_DIR=$(infer_launch_dir "$TRACE_DIR")

if [ -z "$PROCESSED_OUT" ]; then
    PROCESSED_OUT="$TRACE_DIR/processed_sample_workdirs.tsv"
fi
PROCESSED_OUT=$(canonical_output_path "$PROCESSED_OUT")

if [ -n "$WORK_ROOT" ]; then
    [ -d "$WORK_ROOT" ] || fail "--work-root does not exist: $WORK_ROOT"
    WORK_ROOT_CANONICAL=$(canonical_existing_path "$WORK_ROOT")
    WORK_ROOT_EXPLICIT=1
else
    WORK_ROOT_CANONICAL="$LAUNCH_DIR/work"
fi

TMP_DIR=$(mktemp -d "${TMPDIR:-/tmp}/cleanup_processed_sample_workdirs.XXXXXX")
trap 'rm -rf "$TMP_DIR"' EXIT

FINAL_ROWS="$TMP_DIR/final_rows.tsv"
ALL_WORKDIRS="$TMP_DIR/all_workdirs.tsv"
ELIGIBLE_SAMPLES="$TMP_DIR/eligible_samples.tsv"
RAW_WORKDIRS="$TMP_DIR/raw_workdirs.tsv"
SAFE_WORKDIRS="$TMP_DIR/safe_workdirs.tsv"
UNIQUE_PATHS="$TMP_DIR/unique_paths.txt"

parse_trace "$FINAL_ROWS" "$ALL_WORKDIRS"

if [ ! -s "$FINAL_ROWS" ]; then
    log "No sample trace rows matched the requested scope."
    printf 'sra\tsrr\tworkdir_count\tdeleted_default\tdry_run\ttrace_file\n' > "$PROCESSED_OUT"
    exit 0
fi

write_eligible_samples "$FINAL_ROWS" "$ELIGIBLE_SAMPLES"
if [ ! -s "$ELIGIBLE_SAMPLES" ]; then
    log "No processed samples found in trace."
    printf 'sra\tsrr\tworkdir_count\tdeleted_default\tdry_run\ttrace_file\n' > "$PROCESSED_OUT"
    exit 0
fi

filter_eligible_workdirs "$ELIGIBLE_SAMPLES" "$ALL_WORKDIRS" "$RAW_WORKDIRS"
collect_safe_workdirs "$RAW_WORKDIRS" "$SAFE_WORKDIRS"
write_processed_samples "$ELIGIBLE_SAMPLES" "$SAFE_WORKDIRS"

if [ ! -s "$SAFE_WORKDIRS" ]; then
    log "No safe workdirs matched processed samples."
    exit 0
fi

cleanup_workdirs "$SAFE_WORKDIRS" "$UNIQUE_PATHS"
