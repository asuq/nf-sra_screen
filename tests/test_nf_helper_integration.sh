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

    grep -F -- "$pattern" "$file" >/dev/null 2>&1 || fail "expected $file to contain $pattern"
}

assert_path_exists() {
    [ -e "$1" ] || fail "expected path to exist: $1"
}

assert_path_missing() {
    [ ! -e "$1" ] || fail "expected path to be absent: $1"
}

assert_file_omits() {
    local file=$1
    local pattern=$2

    ! grep -F -- "$pattern" "$file" >/dev/null 2>&1 || fail "expected $file to omit $pattern"
}

assert_local_module_run_local() {
    # Assert that a local helper module opts into nf-helper login-node policy.
    local module_file=$1

    assert_file_contains "$module_file" "label 'run_local'"
}

assert_site_local_checks() {
    # Assert that lightweight metadata and Sandpiper checks use a bounded local pool.
    local config_file=$1

    assert_file_contains "$config_file" "withName: 'DOWNLOAD_SRA_METADATA|RESOLVE_SRR_METADATA|SANDPIPER' {"
    assert_file_contains "$config_file" "executor = 'local'"
    assert_file_contains "$config_file" '$local {'
    assert_file_contains "$config_file" "cpus = 2"
    assert_file_contains "$config_file" "memory = 32.GB"
}

assert_process_block_contains() {
    # Assert within one named process block rather than matching another selector.
    local config_file=$1
    local selector=$2
    local expected=$3
    local block

    block="$(
        awk -v selector="$selector" '
        index($0, selector) {
            capture = 1
        }
        capture {
            print
            opening_line = $0
            closing_line = $0
            depth += gsub(/\{/, "{", opening_line)
            depth -= gsub(/\}/, "}", closing_line)
            if (depth == 0) {
                exit
            }
        }
        ' "$config_file"
    )"

    [ -n "$block" ] || fail "expected $config_file to contain selector $selector"
    grep -F -- "$expected" <<<"$block" >/dev/null 2>&1 \
        || fail "expected $selector in $config_file to contain $expected"
}

assert_download_srr_runs_local() {
    local config_file=$1

    assert_process_block_contains "$config_file" "withName: DOWNLOAD_SRR {" "executor = 'local'"
    assert_process_block_contains "$config_file" "withName: DOWNLOAD_SRR {" "scratch = false"
}

assert_file_contains "$REPO_ROOT/.gitmodules" "path = external/nf-helper"
assert_file_contains "$REPO_ROOT/.gitmodules" "url = https://github.com/asuq/nf-helper.git"
assert_path_exists "$REPO_ROOT/external/nf-helper/conf/sites/oist.config"
assert_path_exists "$REPO_ROOT/external/nf-helper/conf/sites/gwdg.config"
assert_path_exists "$REPO_ROOT/external/nf-helper/conf/sites/marmic.config"
assert_path_exists "$REPO_ROOT/external/nf-helper/conf/sites/viper-cpu.config"
assert_path_exists "$REPO_ROOT/external/nf-helper/helpers/cleanup_processed_sample_workdirs.sh"
assert_path_exists "$REPO_ROOT/external/nf-helper/helpers/gwdg_promote_2h_qos.sh"
assert_path_missing "$REPO_ROOT/helpers/cleanup_processed_sra_workdirs.sh"
assert_file_contains "$REPO_ROOT/external/nf-helper/conf/sites/viper-cpu.config" "viper_slurm_queue_size = 250"
assert_file_contains "$REPO_ROOT/external/nf-helper/conf/sites/viper-cpu.config" "queueSize = params.viper_slurm_queue_size as int"

assert_file_contains "$REPO_ROOT/conf/oist.config" "external/nf-helper/conf/sites/oist.config"
assert_file_contains "$REPO_ROOT/conf/gwdg.config" "external/nf-helper/conf/sites/gwdg.config"
assert_file_contains "$REPO_ROOT/conf/marmic.config" "external/nf-helper/conf/sites/marmic.config"
assert_file_contains "$REPO_ROOT/conf/viper-cpu.config" "external/nf-helper/conf/sites/viper-cpu.config"
assert_site_local_checks "$REPO_ROOT/conf/oist.config"
assert_site_local_checks "$REPO_ROOT/conf/gwdg.config"
assert_site_local_checks "$REPO_ROOT/conf/marmic.config"
assert_download_srr_runs_local "$REPO_ROOT/conf/oist.config"
assert_download_srr_runs_local "$REPO_ROOT/conf/gwdg.config"
assert_download_srr_runs_local "$REPO_ROOT/conf/marmic.config"
assert_download_srr_runs_local "$REPO_ROOT/conf/viper-cpu.config"
assert_file_contains "$REPO_ROOT/conf/viper-cpu.config" "'viper-cpu' {"
assert_file_contains "$REPO_ROOT/conf/viper-cpu.config" "withName: VALIDATE_TAXA"
assert_file_contains "$REPO_ROOT/conf/viper-cpu.config" "memory = 16.GB"
assert_file_contains "$REPO_ROOT/conf/viper-cpu.config" "withName: SINGLEM"
assert_file_contains "$REPO_ROOT/conf/viper-cpu.config" "clusterOptions = '--export=ALL --ntasks-per-core=2'"
assert_file_contains "$REPO_ROOT/conf/viper-cpu.config" "cpus = { Math.min(64 * task.attempt, params.max_cpus as int) }"
assert_file_contains "$REPO_ROOT/conf/viper-cpu.config" "def requestedTime = task.attempt == 1 ? 4.h : task.attempt == 2 ? 12.h : 24.h"
assert_file_contains "$REPO_ROOT/conf/viper-cpu.config" "[requestedTime, params.resource_time_limit()].min()"
assert_file_contains "$REPO_ROOT/modules/local/singlem/main.nf" '--cpus ${task.cpus}'
assert_file_contains "$REPO_ROOT/bin/run_singlem.sh" '--threads "$cpus"'
assert_file_contains "$REPO_ROOT/conf/viper-cpu.config" "DOWNLOAD_SRA_METADATA"
assert_file_contains "$REPO_ROOT/conf/viper-cpu.config" "RESOLVE_SRR_METADATA"
assert_file_contains "$REPO_ROOT/conf/viper-cpu.config" "SANDPIPER"
assert_file_contains "$REPO_ROOT/nextflow.config" "download_srr_max_forks = 2"
assert_process_block_contains "$REPO_ROOT/nextflow.config" "withName: DOWNLOAD_SRR {" "maxForks = params.download_srr_max_forks as int"
assert_process_block_contains "$REPO_ROOT/nextflow.config" "withName: DOWNLOAD_SRR {" "cpus   = 1"
assert_process_block_contains "$REPO_ROOT/nextflow.config" "withName: DOWNLOAD_SRR {" "memory = { [ 8.GB, params.resource_memory_limit() ].min() }"
assert_process_block_contains "$REPO_ROOT/nextflow.config" "withName: DOWNLOAD_SRR {" "time   = { params.resource_time_limit() }"
assert_file_contains "$REPO_ROOT/nextflow.config" "includeConfig \"\${projectDir}/conf/oist.config\""
assert_file_contains "$REPO_ROOT/nextflow.config" "includeConfig \"\${projectDir}/conf/gwdg.config\""
assert_file_contains "$REPO_ROOT/nextflow.config" "includeConfig \"\${projectDir}/conf/marmic.config\""
assert_file_contains "$REPO_ROOT/nextflow.config" "includeConfig \"\${projectDir}/conf/viper-cpu.config\""
assert_file_contains "$REPO_ROOT/nextflow.config" "cleanup = false"
assert_file_omits "$REPO_ROOT/nextflow.config" "workDir ="
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
