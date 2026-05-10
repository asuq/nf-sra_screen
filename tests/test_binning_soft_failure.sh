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

    grep -Fq -- "$expected" "$file" || fail "expected '$expected' in $file"
}

assert_file_has_line() {
    # Assert that a file contains an exact line.
    local file=$1
    local expected=$2

    grep -Fxq -- "$expected" "$file" || fail "expected line '$expected' in $file"
}

assert_file_missing() {
    # Assert that a file does not exist.
    [ ! -e "$1" ] || fail "expected file to be absent: $1"
}

append_fasta_record() {
    # Append one synthetic FASTA record of the requested length.
    local file=$1
    local record_id=$2
    local length=$3

    printf '>%s\n' "$record_id" >> "$file"
    awk -v n="$length" 'BEGIN {
      for (i = 0; i < n; i++) {
        printf "A"
      }
      printf "\n"
    }' >> "$file"
}

write_eligible_comebin_assembly() {
    # Write a tiny assembly with enough COMEBin-eligible contigs.
    local assembly=$1

    : > "$assembly"
    append_fasta_record "$assembly" "contig1" 1001
    append_fasta_record "$assembly" "contig2" 1001
}

write_mixed_comebin_assembly() {
    # Write an assembly with eligible and short contigs.
    local assembly=$1

    : > "$assembly"
    append_fasta_record "$assembly" "kept1" 1001
    append_fasta_record "$assembly" "short1" 500
    append_fasta_record "$assembly" "kept2" 1002
}

write_insufficient_comebin_assembly() {
    # Write an assembly with fewer than two COMEBin-eligible contigs.
    local assembly=$1

    : > "$assembly"
    append_fasta_record "$assembly" "kept1" 1001
    append_fasta_record "$assembly" "short1" 500
}

write_fake_comebin() {
    # Write a fake COMEBin executable that always fails.
    local fake_bin=$1
    local exit_code=${2:-42}

    mkdir -p "$fake_bin"
    cat > "$fake_bin/run_comebin.sh" <<'EOF'
#!/usr/bin/env bash
if [[ -n "${FAKE_COMEBIN_ARGS:-}" ]]; then
  printf '%s\n' "$@" > "$FAKE_COMEBIN_ARGS"
fi

output_dir=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    -o) output_dir="$2"; shift 2 ;;
    *) shift ;;
  esac
done

if [[ "${FAKE_COMEBIN_EXIT:-42}" == "0" ]]; then
  mkdir -p "${output_dir}/comebin_res/comebin_res_bins"
  printf '>kept1\nACGT\n' > "${output_dir}/comebin_res/comebin_res_bins/bin.1.fa"
  exit 0
fi

printf 'fake COMEBin failure\n' >&2
exit "${FAKE_COMEBIN_EXIT:-42}"
EOF
    chmod +x "$fake_bin/run_comebin.sh"
    printf '%s\n' "$exit_code" > "$fake_bin/exit_code"
}

run_comebin_in_workdir() {
    # Run run_comebin_nf.sh in an isolated directory with an existing assembly.
    local attempt=$1
    local work_dir=$2
    local fake_bin=$3

    mkdir -p "$work_dir"
    : > "$work_dir/assembly.bam"

    (
        cd "$work_dir"
        FAKE_COMEBIN_EXIT="$(cat "$fake_bin/exit_code")" \
        FAKE_COMEBIN_ARGS="${FAKE_COMEBIN_ARGS:-}" \
        PATH="$fake_bin:$REPO_ROOT/bin:$PATH" \
            run_comebin_nf.sh \
                --assembly assembly.fasta \
                --bam assembly.bam \
                --cpus 1 \
                --attempt "$attempt" \
                --max-retries 1
    )
}

run_comebin_attempt() {
    # Run run_comebin_nf.sh in an isolated directory for one attempt.
    local attempt=$1
    local work_dir=$2
    local fake_bin=$3

    mkdir -p "$work_dir"
    write_eligible_comebin_assembly "$work_dir/assembly.fasta"
    run_comebin_in_workdir "$attempt" "$work_dir" "$fake_bin"
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
    assert_file_contains "$final_work/FAIL.note" "COMEBin: run_comebin.sh failed"
}

test_comebin_timeout_exit_stays_hard_failure() {
    # Scheduler-style timeout exits should not be converted into soft outputs.
    local fake_bin="$TMP_ROOT/fake_timeout_bin"
    local timeout_work="$TMP_ROOT/comebin_timeout"

    write_fake_comebin "$fake_bin" 143

    if run_comebin_attempt 2 "$timeout_work" "$fake_bin" > "$timeout_work.log" 2>&1; then
        fail "COMEBin timeout attempt unexpectedly soft-succeeded"
    fi

    [ ! -d "$timeout_work/comebin" ] || fail "timeout produced a soft COMEBin directory"
    assert_file_empty "$timeout_work/comebin.contig2bin.tsv"
}

test_comebin_batch_uses_filtered_contigs() {
    # COMEBin should receive a batch size based on eligible contigs only.
    local fake_bin="$TMP_ROOT/fake_success_bin"
    local work_dir="$TMP_ROOT/comebin_filtered_batch"
    local args_log="$TMP_ROOT/comebin_args.log"

    write_fake_comebin "$fake_bin" 0
    mkdir -p "$work_dir"
    write_mixed_comebin_assembly "$work_dir/assembly.fasta"

    if ! FAKE_COMEBIN_ARGS="$args_log" run_comebin_in_workdir 1 "$work_dir" "$fake_bin" > "$work_dir.log" 2>&1; then
        sed -n '1,160p' "$work_dir.log" >&2
        fail "COMEBin filtered batch run failed"
    fi

    assert_file_contains "$work_dir.log" "COMEBin: using 2/3 contigs longer than 1000 bp with batch size 2"
    assert_file_has_line "$args_log" "-b"
    assert_file_has_line "$args_log" "2"
    assert_file_empty "$work_dir/FAIL.note"
    assert_file_contains "$work_dir/comebin.contig2bin.tsv" "kept1"
}

test_comebin_skips_insufficient_eligible_contigs() {
    # COMEBin should skip deterministically when too few eligible contigs remain.
    local fake_bin="$TMP_ROOT/fake_skip_bin"
    local work_dir="$TMP_ROOT/comebin_skip"
    local args_log="$TMP_ROOT/comebin_skip_args.log"

    write_fake_comebin "$fake_bin" 0
    mkdir -p "$work_dir"
    write_insufficient_comebin_assembly "$work_dir/assembly.fasta"

    if ! FAKE_COMEBIN_ARGS="$args_log" run_comebin_in_workdir 1 "$work_dir" "$fake_bin" > "$work_dir.log" 2>&1; then
        sed -n '1,160p' "$work_dir.log" >&2
        fail "COMEBin insufficient-contig skip failed"
    fi

    assert_dir_exists "$work_dir/comebin"
    assert_file_empty "$work_dir/comebin.contig2bin.tsv"
    assert_file_contains "$work_dir/FAIL.note" "COMEBin: skipped because only 1 contig(s) are longer than 1000 bp"
    assert_file_missing "$args_log"
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
    test_comebin_timeout_exit_stays_hard_failure
    test_comebin_batch_uses_filtered_contigs
    test_comebin_skips_insufficient_eligible_contigs
    test_partial_binner_group_reaches_refiner_input

    printf 'PASS: binning soft-failure recovery tests\n'
}

main "$@"
