#!/usr/bin/env bash

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." >/dev/null 2>&1 && pwd -P)"
HELPER="$REPO_ROOT/helpers/cleanup_processed_sra_workdirs.sh"
TMP_ROOT=""

fail() {
    # Print an assertion failure and exit.
    printf 'FAIL: %s\n' "$*" >&2
    exit 1
}

assert_dir_exists() {
    # Assert that a directory exists.
    [ -d "$1" ] || fail "expected directory to exist: $1"
}

assert_dir_missing() {
    # Assert that a directory does not exist.
    [ ! -e "$1" ] || fail "expected path to be removed: $1"
}

assert_file_contains() {
    # Assert that a file contains a fixed string.
    local file=$1
    local pattern=$2

    grep -F "$pattern" "$file" >/dev/null 2>&1 || \
        fail "expected '$file' to contain '$pattern'"
}

setup_fixture() {
    # Create a run directory with processed samples, trace rows, and workdirs.
    local run_dir=$1

    mkdir -p \
        "$run_dir/execution-reports" \
        "$run_dir/work/aa/bb" \
        "$run_dir/work/cc/dd" \
        "$run_dir/work/ee/ff" \
        "$run_dir/work/sraonly" \
        "$run_dir/outside"

    printf 'a\n' > "$run_dir/work/aa/bb/file.txt"
    printf 'b\n' > "$run_dir/work/cc/dd/file.txt"
    printf 'c\n' > "$run_dir/work/ee/ff/file.txt"
    printf 'd\n' > "$run_dir/work/sraonly/file.txt"
    printf 'e\n' > "$run_dir/outside/file.txt"

    cat > "$run_dir/.processed_summary.tsv" <<EOF
SRA1	SRR1
SRA2	SRR2
EOF

    cat > "$run_dir/execution-reports/trace.tsv" <<EOF
process	status	tag	workdir
DOWNLOAD	COMPLETED	SRA1:SRR1	$run_dir/work/aa/bb
ASSEMBLE	COMPLETED	SRA1:SRR1:metaspades	$run_dir/work/aa/bb
MAP	COMPLETED	SRA2:SRR2	work/cc/dd
OTHER	COMPLETED	SRA3:SRR3	$run_dir/work/ee/ff
METADATA	COMPLETED	SRA1	$run_dir/work/sraonly
BAD	COMPLETED	SRA1:SRR1	$run_dir/outside
MISSING	COMPLETED	SRA2:SRR2	$run_dir/work/missing
EOF
}

test_dry_run_keeps_workdirs() {
    # Verify dry-run mode reports matches without deleting anything.
    local run_dir="$TMP_ROOT/dry_run"
    local log_file="$TMP_ROOT/dry_run.log"

    setup_fixture "$run_dir"

    "$HELPER" "$run_dir" --dry-run > "$log_file" 2>&1

    assert_dir_exists "$run_dir/work/aa/bb"
    assert_dir_exists "$run_dir/work/cc/dd"
    assert_dir_exists "$run_dir/work/ee/ff"
    assert_dir_exists "$run_dir/work/sraonly"
    assert_dir_exists "$run_dir/outside"
    assert_file_contains "$log_file" "DRY-RUN summary: matched_workdirs=2"
}

test_default_deletes_only_safe_processed_workdirs() {
    # Verify default mode deletes only matching workdirs under RUN_DIR/work.
    local run_dir="$TMP_ROOT/delete"
    local log_file="$TMP_ROOT/delete.log"

    setup_fixture "$run_dir"

    "$HELPER" "$run_dir" > "$log_file" 2>&1

    assert_dir_missing "$run_dir/work/aa/bb"
    assert_dir_missing "$run_dir/work/cc/dd"
    assert_dir_exists "$run_dir/work/ee/ff"
    assert_dir_exists "$run_dir/work/sraonly"
    assert_dir_exists "$run_dir/outside"
    assert_file_contains "$log_file" "Skipping suspicious workdir outside RUN_DIR/work"
    assert_file_contains "$log_file" "Skipping missing workdir from trace"
    assert_file_contains "$log_file" "Cleanup summary: deleted_workdirs=2"
}

test_missing_trace_header_fails() {
    # Verify malformed trace files fail clearly.
    local run_dir="$TMP_ROOT/bad_trace"
    local log_file="$TMP_ROOT/bad_trace.log"

    mkdir -p "$run_dir/work" "$run_dir/execution-reports"
    printf 'SRA1\tSRR1\n' > "$run_dir/.processed_summary.tsv"
    cat > "$run_dir/execution-reports/trace.tsv" <<EOF
process	tag
DOWNLOAD	SRA1:SRR1
EOF

    if "$HELPER" "$run_dir" > "$log_file" 2>&1; then
        fail "malformed trace command unexpectedly succeeded"
    fi

    assert_file_contains "$log_file" "trace TSV must include tag and workdir headers"
}

main() {
    # Run all fixture tests.
    TMP_ROOT=$(mktemp -d "${TMPDIR:-/tmp}/cleanup_processed_sra_workdirs_test.XXXXXX")
    trap 'rm -rf "$TMP_ROOT"' EXIT

    test_dry_run_keeps_workdirs
    test_default_deletes_only_safe_processed_workdirs
    test_missing_trace_header_fails

    printf 'cleanup_processed_sra_workdirs tests passed\n'
}

main "$@"
