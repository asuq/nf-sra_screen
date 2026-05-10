#!/usr/bin/env bash

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." >/dev/null 2>&1 && pwd -P)"
TMP_ROOT=""

fail() {
    # Print an assertion failure and exit.
    printf 'FAIL: %s\n' "$*" >&2
    exit 1
}

assert_file_exists() {
    # Assert that a file exists.
    [ -f "$1" ] || fail "expected file to exist: $1"
}

assert_dir_exists() {
    # Assert that a directory exists.
    [ -d "$1" ] || fail "expected directory to exist: $1"
}

assert_file_empty() {
    # Assert that a file exists and has zero bytes.
    assert_file_exists "$1"
    [ ! -s "$1" ] || fail "expected file to be empty: $1"
}

assert_file_contains() {
    # Assert that a file contains a fixed string.
    local file=$1
    local expected=$2

    grep -Fq "$expected" "$file" || fail "expected '$expected' in $file"
}

write_fake_comebin() {
    # Write a fake COMEBin executable that always fails.
    local fake_bin=$1

    mkdir -p "$fake_bin"
    cat > "$fake_bin/run_comebin.sh" <<'EOF'
#!/usr/bin/env bash
printf 'fake COMEBin failure\n' >&2
exit 42
EOF
    chmod +x "$fake_bin/run_comebin.sh"
}

run_comebin_attempt() {
    # Run run_comebin_nf.sh in an isolated directory for one attempt.
    local attempt=$1
    local work_dir=$2
    local fake_bin=$3

    mkdir -p "$work_dir"
    printf '>contig1\nACGT\n' > "$work_dir/assembly.fasta"
    : > "$work_dir/assembly.bam"

    (
        cd "$work_dir"
        PATH="$fake_bin:$REPO_ROOT/bin:$PATH" \
            run_comebin_nf.sh \
                --assembly assembly.fasta \
                --bam assembly.bam \
                --cpus 1 \
                --attempt "$attempt" \
                --max-retries 1
    )
}

test_comebin_retries_before_soft_failure() {
    # COMEBin should fail before the final retry and soft-succeed afterwards.
    local fake_bin="$TMP_ROOT/fake_bin"
    local retry_work="$TMP_ROOT/comebin_retry"
    local final_work="$TMP_ROOT/comebin_final"

    write_fake_comebin "$fake_bin"

    if run_comebin_attempt 1 "$retry_work" "$fake_bin" > "$retry_work.log" 2>&1; then
        fail "COMEBin retry attempt unexpectedly succeeded"
    fi
    assert_file_contains "$retry_work.log" "COMEBin: run_comebin.sh failed"

    if ! run_comebin_attempt 2 "$final_work" "$fake_bin" > "$final_work.log" 2>&1; then
        sed -n '1,120p' "$final_work.log" >&2
        fail "COMEBin final soft-failure attempt failed"
    fi

    assert_dir_exists "$final_work/comebin"
    assert_file_empty "$final_work/comebin.contig2bin.tsv"
    assert_file_contains "$final_work/comebin.note" "COMEBin: run_comebin.sh failed"
}

write_partial_group_fixture() {
    # Write a tiny Nextflow fixture for one missing binner result.
    local fixture_file=$1

    cat > "$fixture_file" <<'EOF'
nextflow.enable.dsl=2

workflow {
  def metabat_map = file(params.metabat_map, checkIfExists: true)
  def metabat_dir = file(params.metabat_dir, checkIfExists: true)
  def metabat_note = file(params.metabat_note, checkIfExists: true)

  def binner_plan_by_sample = channel.of(
    tuple(['SRA1', 'SRR1', 'short', 'metaspades'], 2)
  )

  def binner_results = channel.of(
    tuple(
      'SRA1',
      'SRR1',
      'ILLUMINA',
      'NovaSeq',
      'WGS',
      'short',
      'metaspades',
      'metabat',
      metabat_dir,
      metabat_map,
      metabat_note
    )
  )

  binner_results
    .map { sra, srr, platform, model, strategy, read_type, assembler, tool, bin_dir, contig2bin, note_path ->
      tuple([sra, srr, read_type, assembler], [tool, contig2bin])
    }
    .combine(binner_plan_by_sample, by: 0)
    .map { sample_join_key, binner_map_payload, expected_count ->
      tuple(groupKey(sample_join_key, expected_count as int), binner_map_payload)
    }
    .groupTuple(remainder: true)
    .map { grouped_sample_key, binner_map_entries ->
      tuple(grouped_sample_key.target, binner_map_entries)
    }
    .view { sample_join_key, binner_map_entries ->
      def tools = binner_map_entries.collect { entry -> entry[0] }.join(',')
      def maps = binner_map_entries.collect { entry -> entry[1].getName() }.join(',')
      "GROUP ${sample_join_key.join(':')} count=${binner_map_entries.size()} tools=${tools} maps=${maps}"
    }
}
EOF
}

test_partial_binner_group_reaches_refiner_input() {
    # A missing COMEBin result should not prevent a partial refiner input.
    local fixture_file="$TMP_ROOT/partial_group_fixture.nf"
    local workdir="$TMP_ROOT/partial_group_work"
    local log_file="$TMP_ROOT/partial_group.log"
    local empty_config="$TMP_ROOT/empty_nextflow.config"
    local metabat_dir="$TMP_ROOT/metabat"
    local metabat_map="$TMP_ROOT/metabat.contig2bin.tsv"
    local metabat_note="$TMP_ROOT/metabat.note"

    mkdir -p "$metabat_dir"
    printf 'contig1\tmetabatbin.1.fa\n' > "$metabat_map"
    : > "$metabat_note"
    : > "$empty_config"
    write_partial_group_fixture "$fixture_file"

    if ! NXF_OFFLINE=true nextflow -C "$empty_config" run "$fixture_file" \
        -ansi-log false \
        -work-dir "$workdir" \
        --metabat_dir "$metabat_dir" \
        --metabat_map "$metabat_map" \
        --metabat_note "$metabat_note" \
        > "$log_file" 2>&1; then
        sed -n '1,160p' "$log_file" >&2
        fail "partial binner group fixture failed"
    fi

    assert_file_contains "$log_file" "GROUP SRA1:SRR1:short:metaspades count=1"
    assert_file_contains "$log_file" "tools=metabat"
    assert_file_contains "$log_file" "maps=metabat.contig2bin.tsv"
}

main() {
    # Run all soft-failure recovery tests.
    TMP_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/nf-sra-binning-soft-failure.XXXXXX")"
    trap 'rm -rf "$TMP_ROOT"' EXIT

    test_comebin_retries_before_soft_failure
    test_partial_binner_group_reaches_refiner_input

    printf 'PASS: binning soft-failure recovery tests\n'
}

main "$@"
