#!/usr/bin/env bash
#
# Delete Nextflow work directories for samples already processed by
# watch_and_transfer.sh.
#
# The helper reads RUN_DIR/.processed_summary.tsv and a Nextflow trace file,
# matches processed (sra, srr) pairs to trace tags, and deletes only matching
# work directories under RUN_DIR/work.

set -euo pipefail

DRY_RUN=0
RUN_DIR=""
STATE_FILE=""
TRACE_FILE=""

usage() {
    # Print command-line usage.
    cat <<EOF
Usage: $(basename "$0") RUN_DIR [options]

Options:
  --trace TRACE_TSV       Trace file to read (default: RUN_DIR/execution-reports/trace.tsv, then RUN_DIR/trace.tsv)
  --state-file STATE_TSV  Processed sample TSV (default: RUN_DIR/.processed_summary.tsv)
  --dry-run               Print matched work directories and sizes without deleting
  -h, --help              Show this help message

The trace file must be tab-separated and include 'tag' and 'workdir' headers.
Only tags shaped like SRA:SRR or SRA:SRR:* are eligible for cleanup.
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

default_trace_file() {
    # Print the default trace path for RUN_DIR.
    local run_dir=$1
    local candidate

    candidate="$run_dir/execution-reports/trace.tsv"
    if [ -r "$candidate" ]; then
        printf '%s\n' "$candidate"
        return 0
    fi

    candidate="$run_dir/trace.tsv"
    if [ -r "$candidate" ]; then
        printf '%s\n' "$candidate"
        return 0
    fi

    return 1
}

normalise_trace_workdir() {
    # Print an absolute candidate workdir path from a trace workdir field.
    local workdir=$1

    case "$workdir" in
        /*)
            printf '%s\n' "$workdir"
            ;;
        *)
            printf '%s/%s\n' "$RUN_DIR" "$workdir"
            ;;
    esac
}

extract_processed_pairs() {
    # Extract unique non-empty (sra, srr) pairs from the state file.
    local output=$1

    awk -F'\t' '
        $1 != "" && $2 != "" { print $1 "\t" $2 }
    ' "$STATE_FILE" | sort -u > "$output"
}

extract_trace_workdirs() {
    # Extract trace workdir fields whose tags match processed pairs.
    local pairs_file=$1
    local output=$2

    awk -F'\t' -v pairs_file="$pairs_file" '
        BEGIN {
            while ((getline line < pairs_file) > 0) {
                split(line, fields, FS)
                if (fields[1] != "" && fields[2] != "") {
                    processed[fields[1] FS fields[2]] = 1
                }
            }
            close(pairs_file)
        }
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                if ($i == "tag") {
                    tag_col = i
                }
                if ($i == "workdir") {
                    workdir_col = i
                }
            }
            if (!tag_col || !workdir_col) {
                print "ERROR: trace TSV must include tag and workdir headers" > "/dev/stderr"
                exit 2
            }
            next
        }
        {
            tag = $tag_col
            workdir = $workdir_col
            if (tag == "" || workdir == "") {
                next
            }

            part_count = split(tag, parts, ":")
            if (part_count < 2) {
                next
            }

            key = parts[1] FS parts[2]
            if (key in processed) {
                print workdir
            }
        }
        END {
            if (NR == 0) {
                print "ERROR: trace TSV is empty" > "/dev/stderr"
                exit 2
            }
        }
    ' "$TRACE_FILE" > "$output"
}

collect_safe_workdirs() {
    # Resolve, fence, and deduplicate matching work directories.
    local raw_paths=$1
    local output=$2
    local path
    local candidate
    local resolved

    : > "$output"

    while IFS= read -r path; do
        [ -n "$path" ] || continue

        candidate=$(normalise_trace_workdir "$path")
        if ! resolved=$(canonical_existing_path "$candidate"); then
            log "Skipping missing workdir from trace: $candidate"
            continue
        fi

        case "$resolved" in
            "$WORK_DIR_CANONICAL"/*)
                printf '%s\n' "$resolved" >> "$output"
                ;;
            *)
                log "Skipping suspicious workdir outside RUN_DIR/work: $resolved"
                ;;
        esac
    done < "$raw_paths"

    sort -u "$output" -o "$output"
}

directory_size_kib() {
    # Print a directory size estimate in KiB.
    du -sk "$1" 2>/dev/null | awk '{ print $1 }'
}

cleanup_workdirs() {
    # Delete or report the fenced work directories.
    local safe_paths=$1
    local count=0
    local total_kib=0
    local size_kib
    local path

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
    done < "$safe_paths"

    if [ "$DRY_RUN" -eq 1 ]; then
        log "DRY-RUN summary: matched_workdirs=$count total_size=${total_kib} KiB"
    else
        log "Cleanup summary: deleted_workdirs=$count total_size=${total_kib} KiB"
    fi
}

while [ "$#" -gt 0 ]; do
    case "$1" in
        --trace)
            [ "$#" -ge 2 ] || fail "--trace requires a value"
            TRACE_FILE=$2
            shift 2
            ;;
        --trace=*)
            TRACE_FILE=${1#*=}
            shift
            ;;
        --state-file)
            [ "$#" -ge 2 ] || fail "--state-file requires a value"
            STATE_FILE=$2
            shift 2
            ;;
        --state-file=*)
            STATE_FILE=${1#*=}
            shift
            ;;
        --dry-run)
            DRY_RUN=1
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
            if [ -n "$RUN_DIR" ]; then
                fail "unexpected positional argument: $1"
            fi
            RUN_DIR=$1
            shift
            ;;
    esac
done

[ -n "$RUN_DIR" ] || { usage >&2; exit 1; }
[ -d "$RUN_DIR" ] || fail "RUN_DIR does not exist: $RUN_DIR"

require_command awk
require_command basename
require_command dirname
require_command du
require_command mktemp
require_command rm
require_command sort

RUN_DIR=$(canonical_existing_path "$RUN_DIR")
WORK_DIR="$RUN_DIR/work"

if [ -z "$STATE_FILE" ]; then
    STATE_FILE="$RUN_DIR/.processed_summary.tsv"
fi

if [ -z "$TRACE_FILE" ]; then
    if ! TRACE_FILE=$(default_trace_file "$RUN_DIR"); then
        fail "trace file not found; expected $RUN_DIR/execution-reports/trace.tsv or $RUN_DIR/trace.tsv"
    fi
fi

[ -r "$STATE_FILE" ] || fail "state file is not readable: $STATE_FILE"
[ -r "$TRACE_FILE" ] || fail "trace file is not readable: $TRACE_FILE"

if [ ! -d "$WORK_DIR" ]; then
    log "No RUN_DIR/work directory exists; nothing to clean: $WORK_DIR"
    exit 0
fi

WORK_DIR_CANONICAL=$(canonical_existing_path "$WORK_DIR")

TMP_DIR=$(mktemp -d "${TMPDIR:-/tmp}/cleanup_processed_sra_workdirs.XXXXXX")
trap 'rm -rf "$TMP_DIR"' EXIT

PAIRS_FILE="$TMP_DIR/processed_pairs.tsv"
RAW_WORKDIRS="$TMP_DIR/raw_workdirs.txt"
SAFE_WORKDIRS="$TMP_DIR/safe_workdirs.txt"

extract_processed_pairs "$PAIRS_FILE"
if [ ! -s "$PAIRS_FILE" ]; then
    log "No processed samples found in state file: $STATE_FILE"
    exit 0
fi

extract_trace_workdirs "$PAIRS_FILE" "$RAW_WORKDIRS"
if [ ! -s "$RAW_WORKDIRS" ]; then
    log "No trace workdirs matched processed samples."
    exit 0
fi

collect_safe_workdirs "$RAW_WORKDIRS" "$SAFE_WORKDIRS"
if [ ! -s "$SAFE_WORKDIRS" ]; then
    log "No safe workdirs matched processed samples."
    exit 0
fi

cleanup_workdirs "$SAFE_WORKDIRS"
