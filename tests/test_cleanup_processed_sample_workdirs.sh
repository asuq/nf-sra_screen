#!/usr/bin/env bash

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." >/dev/null 2>&1 && pwd -P)"
HELPER="$REPO_ROOT/helpers/cleanup_processed_sample_workdirs.sh"
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

assert_file_not_contains() {
    # Assert that a file does not contain a fixed string.
    local file=$1
    local pattern=$2

    if grep -F "$pattern" "$file" >/dev/null 2>&1; then
        fail "expected '$file' not to contain '$pattern'"
    fi
}

canonical_existing_path() {
    # Resolve an existing path to an absolute physical path.
    local path=$1
    local dir
    local base

    dir=$(dirname "$path")
    base=$(basename "$path")
    (
        cd "$dir" >/dev/null 2>&1
        printf '%s/%s\n' "$(pwd -P)" "$base"
    )
}

setup_fixture() {
    # Create a run directory with trace rows and workdirs.
    local run_dir=$1

    mkdir -p \
        "$run_dir/execution-reports" \
        "$run_dir/work/sra1_download_fail" \
        "$run_dir/work/sra1_download_ok" \
        "$run_dir/work/sra1_singlem" \
        "$run_dir/work/sra2_download" \
        "$run_dir/work/sra2_singlem" \
        "$run_dir/work/sra3_running" \
        "$run_dir/work/sra4_dup" \
        "$run_dir/work/sra5_sraonly" \
        "$run_dir/work/sra6_retry_fail" \
        "$run_dir/work/sra6_retry_ok" \
        "$run_dir/outside"

    for dir in "$run_dir"/work/* "$run_dir/outside"; do
        printf 'x\n' > "$dir/file.txt"
    done

    cat > "$run_dir/execution-reports/trace.tsv" <<EOF
task_id	process	name	tag	status	hash	attempt	workdir
1	DOWNLOAD	DOWNLOAD (SRA1:SRR1)	SRA1:SRR1	FAILED	hash_a	1	$run_dir/work/sra1_download_fail
2	DOWNLOAD	DOWNLOAD (SRA1:SRR1)	SRA1:SRR1	COMPLETED	hash_a	2	$run_dir/work/sra1_download_ok
3	SINGLEM	SINGLEM (SRA1:SRR1)	SRA1:SRR1	CACHED	hash_b	1	work/sra1_singlem
4	OUTSIDE	OUTSIDE (SRA1:SRR1)	SRA1:SRR1	COMPLETED	hash_c	1	$run_dir/outside
5	MISSING	MISSING (SRA1:SRR1)	SRA1:SRR1	COMPLETED	hash_d	1	$run_dir/work/missing
6	DOWNLOAD	DOWNLOAD (SRA2:SRR2)	SRA2:SRR2	COMPLETED	hash_e	1	$run_dir/work/sra2_download
7	SINGLEM	SINGLEM (SRA2:SRR2)	SRA2:SRR2	FAILED	hash_f	1	$run_dir/work/sra2_singlem
8	SCREEN	SCREEN (SRA3:SRR3)	SRA3:SRR3	RUNNING	hash_g	1	$run_dir/work/sra3_running
9	DUP	DUP (SRA4:SRR4)	SRA4:SRR4	COMPLETED	hash_h	1	$run_dir/work/sra4_dup
10	DUP2	DUP2 (SRA4:SRR4)	SRA4:SRR4	COMPLETED	hash_i	1	$run_dir/work/sra4_dup
11	METADATA	METADATA (SRA5)	SRA5	COMPLETED	hash_j	1	$run_dir/work/sra5_sraonly
12	RETRY	RETRY (SRA6:SRR6)	SRA6:SRR6	FAILED	hash_k	1	$run_dir/work/sra6_retry_fail
13	RETRY	RETRY (SRA6:SRR6)	SRA6:SRR6	COMPLETED	hash_k	2	$run_dir/work/sra6_retry_ok
EOF
}

test_default_deletes_only_complete_samples() {
    # Verify default deletion only for samples with all final rows complete.
    local run_dir="$TMP_ROOT/default"
    local trace_file="$run_dir/execution-reports/trace.tsv"
    local processed_file="$run_dir/execution-reports/processed_sample_workdirs.tsv"
    local log_file="$TMP_ROOT/default.log"
    local trace_file_canonical

    setup_fixture "$run_dir"
    trace_file_canonical=$(canonical_existing_path "$trace_file")

    "$HELPER" "$trace_file" --work-root "$run_dir/work" > "$log_file" 2>&1

    assert_dir_missing "$run_dir/work/sra1_download_fail"
    assert_dir_missing "$run_dir/work/sra1_download_ok"
    assert_dir_missing "$run_dir/work/sra1_singlem"
    assert_dir_exists "$run_dir/work/sra2_download"
    assert_dir_exists "$run_dir/work/sra2_singlem"
    assert_dir_exists "$run_dir/work/sra3_running"
    assert_dir_missing "$run_dir/work/sra4_dup"
    assert_dir_exists "$run_dir/work/sra5_sraonly"
    assert_dir_missing "$run_dir/work/sra6_retry_fail"
    assert_dir_missing "$run_dir/work/sra6_retry_ok"
    assert_dir_exists "$run_dir/outside"

    assert_file_contains "$processed_file" "sra	srr	workdir_count	deleted_default	dry_run	trace_file"
    assert_file_contains "$processed_file" "SRA1	SRR1	3	true	false	$trace_file_canonical"
    assert_file_contains "$processed_file" "SRA4	SRR4	1	true	false	$trace_file_canonical"
    assert_file_contains "$processed_file" "SRA6	SRR6	2	true	false	$trace_file_canonical"
    assert_file_not_contains "$processed_file" "SRA2	SRR2"
    assert_file_not_contains "$processed_file" "SRA3	SRR3"
    assert_file_not_contains "$processed_file" "SRA5"

    assert_file_contains "$log_file" "Skipping sample SRA2	SRR2: final status includes FAILED"
    assert_file_contains "$log_file" "Skipping sample SRA3	SRR3: final status includes RUNNING"
    assert_file_contains "$log_file" "Skipping workdir outside --work-root"
    assert_file_contains "$log_file" "Skipping missing workdir from trace"
    assert_file_contains "$log_file" "Cleanup summary: processed_samples=3 deleted_workdirs=6"
}

test_dry_run_and_filters() {
    # Verify dry-run preserves directories and filters processed output.
    local run_dir="$TMP_ROOT/dry_run"
    local trace_file="$run_dir/execution-reports/trace.tsv"
    local processed_file="$TMP_ROOT/filtered.tsv"
    local log_file="$TMP_ROOT/dry_run.log"
    local trace_file_canonical

    setup_fixture "$run_dir"
    trace_file_canonical=$(canonical_existing_path "$trace_file")

    "$HELPER" "$trace_file" \
        --work-root "$run_dir/work" \
        --sra SRA1 \
        --srr SRR1 \
        --processed-out "$processed_file" \
        --dry-run > "$log_file" 2>&1

    assert_dir_exists "$run_dir/work/sra1_download_fail"
    assert_dir_exists "$run_dir/work/sra1_download_ok"
    assert_dir_exists "$run_dir/work/sra1_singlem"
    assert_file_contains "$processed_file" "SRA1	SRR1	3	false	true	$trace_file_canonical"
    assert_file_not_contains "$processed_file" "SRA4	SRR4"
    assert_file_contains "$log_file" "DRY-RUN would delete"
    assert_file_contains "$log_file" "DRY-RUN summary: processed_samples=1 matched_workdirs=3"
}

test_inferred_work_root_fence() {
    # Verify the default fence is the launch directory work root.
    local run_dir="$TMP_ROOT/inferred"
    local foreign_dir="$TMP_ROOT/foreign/work/sra7_unsafe"
    local trace_file="$run_dir/execution-reports/trace.tsv"
    local processed_file="$run_dir/execution-reports/processed_sample_workdirs.tsv"
    local log_file="$TMP_ROOT/inferred.log"
    local trace_file_canonical

    mkdir -p "$run_dir/execution-reports" "$run_dir/work/sra7_safe" "$foreign_dir"
    printf 'x\n' > "$run_dir/work/sra7_safe/file.txt"
    printf 'x\n' > "$foreign_dir/file.txt"
    cat > "$trace_file" <<EOF
task_id	process	name	tag	status	hash	attempt	workdir
1	SAFE	SAFE (SRA7:SRR7)	SRA7:SRR7	COMPLETED	hash_l	1	$run_dir/work/sra7_safe
2	UNSAFE	UNSAFE (SRA7:SRR7)	SRA7:SRR7	COMPLETED	hash_m	1	$foreign_dir
EOF
    trace_file_canonical=$(canonical_existing_path "$trace_file")

    "$HELPER" "$trace_file" > "$log_file" 2>&1

    assert_dir_missing "$run_dir/work/sra7_safe"
    assert_dir_exists "$foreign_dir"
    assert_file_contains "$processed_file" "SRA7	SRR7	1	true	false	$trace_file_canonical"
    assert_file_contains "$log_file" "Skipping workdir outside inferred work root"
    assert_file_contains "$log_file" "Cleanup summary: processed_samples=1 deleted_workdirs=1"
}

test_missing_trace_header_fails() {
    # Verify malformed trace files fail clearly.
    local run_dir="$TMP_ROOT/bad_trace"
    local trace_file="$run_dir/execution-reports/trace.tsv"
    local log_file="$TMP_ROOT/bad_trace.log"

    mkdir -p "$run_dir/execution-reports"
    cat > "$trace_file" <<EOF
process	tag	status
DOWNLOAD	SRA1:SRR1	COMPLETED
EOF

    if "$HELPER" "$trace_file" > "$log_file" 2>&1; then
        fail "malformed trace command unexpectedly succeeded"
    fi

    assert_file_contains "$log_file" "trace TSV must include tag, status, and workdir headers"
}

main() {
    # Run all fixture tests.
    TMP_ROOT=$(mktemp -d "${TMPDIR:-/tmp}/cleanup_processed_sample_workdirs_test.XXXXXX")
    trap 'rm -rf "$TMP_ROOT"' EXIT

    test_default_deletes_only_complete_samples
    test_dry_run_and_filters
    test_inferred_work_root_fence
    test_missing_trace_header_fails

    printf 'cleanup_processed_sample_workdirs tests passed\n'
}

main "$@"
