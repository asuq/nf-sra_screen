#!/usr/bin/env bash

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." >/dev/null 2>&1 && pwd -P)"

fail() {
    printf 'FAIL: %s\n' "$*" >&2
    exit 1
}

assert_file_contains() {
    local file=$1
    local pattern=$2

    grep -F "$pattern" "$file" >/dev/null 2>&1 || fail "expected $file to contain $pattern"
}

assert_path_exists() {
    [ -e "$1" ] || fail "expected path to exist: $1"
}

assert_path_missing() {
    [ ! -e "$1" ] || fail "expected path to be absent: $1"
}

assert_local_module_run_local() {
    # Assert that a local helper module opts into nf-helper login-node policy.
    local module_file=$1

    assert_file_contains "$module_file" "label 'run_local'"
}

assert_file_contains "$REPO_ROOT/.gitmodules" "path = external/nf-helper"
assert_file_contains "$REPO_ROOT/.gitmodules" "url = https://github.com/asuq/nf-helper.git"
assert_path_exists "$REPO_ROOT/external/nf-helper/conf/sites/oist.config"
assert_path_exists "$REPO_ROOT/external/nf-helper/conf/sites/gwdg.config"
assert_path_exists "$REPO_ROOT/external/nf-helper/conf/sites/marmic.config"
assert_path_exists "$REPO_ROOT/external/nf-helper/helpers/cleanup_processed_sample_workdirs.sh"
assert_path_exists "$REPO_ROOT/external/nf-helper/helpers/gwdg_promote_2h_qos.sh"
assert_path_missing "$REPO_ROOT/helpers/cleanup_processed_sra_workdirs.sh"

assert_file_contains "$REPO_ROOT/conf/oist.config" "external/nf-helper/conf/sites/oist.config"
assert_file_contains "$REPO_ROOT/conf/gwdg.config" "external/nf-helper/conf/sites/gwdg.config"
assert_file_contains "$REPO_ROOT/conf/marmic.config" "external/nf-helper/conf/sites/marmic.config"
assert_file_contains "$REPO_ROOT/nextflow.config" "includeConfig \"\${projectDir}/conf/oist.config\""
assert_file_contains "$REPO_ROOT/nextflow.config" "includeConfig \"\${projectDir}/conf/gwdg.config\""
assert_file_contains "$REPO_ROOT/nextflow.config" "includeConfig \"\${projectDir}/conf/marmic.config\""
assert_local_module_run_local "$REPO_ROOT/modules/local/create_empty_summary/main.nf"
assert_local_module_run_local "$REPO_ROOT/modules/local/create_empty_failure_summary/main.nf"
assert_local_module_run_local "$REPO_ROOT/modules/local/create_assembler_selection_note/main.nf"
assert_local_module_run_local "$REPO_ROOT/modules/local/append_summary/main.nf"
assert_file_contains "$REPO_ROOT/modules/local/dastool/main.nf" "label 'binning'"
assert_file_contains "$REPO_ROOT/modules/local/binette/main.nf" "label 'binning'"
assert_file_contains "$REPO_ROOT/modules/local/comebin_gpu/main.nf" "label 'gpu'"
assert_file_contains "$REPO_ROOT/modules/local/comebin_gpu/main.nf" "label 'binning'"
assert_file_contains "$REPO_ROOT/modules/local/vamb_gpu/main.nf" "label 'gpu'"
assert_file_contains "$REPO_ROOT/modules/local/vamb_gpu/main.nf" "label 'binning'"
assert_file_contains "$REPO_ROOT/modules/local/lorbin_gpu/main.nf" "label 'gpu'"
assert_file_contains "$REPO_ROOT/modules/local/lorbin_gpu/main.nf" "label 'binning'"

bash -n "$REPO_ROOT/helpers/cleanup_processed_sample_workdirs.sh"
bash -n "$REPO_ROOT/helpers/gwdg_promote_2h_qos.sh"
"$REPO_ROOT/helpers/cleanup_processed_sample_workdirs.sh" --help >/dev/null
"$REPO_ROOT/helpers/gwdg_promote_2h_qos.sh" --help >/dev/null

printf 'nf-helper integration tests passed\n'
