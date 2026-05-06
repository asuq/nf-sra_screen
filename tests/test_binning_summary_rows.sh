#!/usr/bin/env bash

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." >/dev/null 2>&1 && pwd -P)"
TMP_ROOT=""
SUMMARY_CSV=""

fail() {
    # Print an assertion failure and exit.
    printf 'FAIL: %s\n' "$*" >&2
    exit 1
}

assert_file_exists() {
    # Assert that a file exists.
    [ -f "$1" ] || fail "expected file to exist: $1"
}

assert_summary_header() {
    # Assert that summary.tsv keeps the expected column order.
    local summary_file=$1
    local expected
    local actual

    expected=$'sra\tsrr\tplatform\tmodel\tstrategy\tread_type\tassembler\tcounts\tnote'
    actual="$(sed -n '1p' "$summary_file")"
    [ "$actual" = "$expected" ] || fail "unexpected header in $summary_file: $actual"
}

assert_data_row_count() {
    # Assert the number of data rows below the header.
    local summary_file=$1
    local expected=$2
    local actual

    actual="$(awk 'END { print (NR > 0 ? NR - 1 : 0) }' "$summary_file")"
    [ "$actual" = "$expected" ] || fail "expected $expected data rows in $summary_file, saw $actual"
}

assert_key_row_count() {
    # Assert row cardinality for one summary key.
    local summary_file=$1
    local sra=$2
    local srr=$3
    local read_type=$4
    local assembler=$5
    local expected=$6
    local actual

    actual="$(
        awk -F '\t' \
            -v sra="$sra" \
            -v srr="$srr" \
            -v read_type="$read_type" \
            -v assembler="$assembler" \
            'NR > 1 && $1 == sra && $2 == srr && $6 == read_type && $7 == assembler { count++ }
             END { print count + 0 }' \
            "$summary_file"
    )"
    [ "$actual" = "$expected" ] || \
        fail "expected $expected rows for $sra:$srr:$read_type:$assembler in $summary_file, saw $actual"
}

summary_note_for_key() {
    # Print the note field for one summary key.
    local summary_file=$1
    local sra=$2
    local srr=$3
    local read_type=$4
    local assembler=$5

    awk -F '\t' \
        -v sra="$sra" \
        -v srr="$srr" \
        -v read_type="$read_type" \
        -v assembler="$assembler" \
        'NR > 1 && $1 == sra && $2 == srr && $6 == read_type && $7 == assembler {
           print $9
           found = 1
           exit
         }
         END {
           if (!found) {
             exit 2
           }
         }' \
        "$summary_file"
}

assert_key_note_blank() {
    # Assert that one summary key has a blank note field.
    local summary_file=$1
    local sra=$2
    local srr=$3
    local read_type=$4
    local assembler=$5
    local note

    note="$(summary_note_for_key "$summary_file" "$sra" "$srr" "$read_type" "$assembler")" || \
        fail "missing row for $sra:$srr:$read_type:$assembler in $summary_file"
    [ -z "$note" ] || fail "expected blank note for $sra:$srr:$read_type:$assembler, saw: $note"
}

assert_key_note_contains() {
    # Assert that one summary key's note contains a fixed string.
    local summary_file=$1
    local sra=$2
    local srr=$3
    local read_type=$4
    local assembler=$5
    local expected=$6
    local note

    note="$(summary_note_for_key "$summary_file" "$sra" "$srr" "$read_type" "$assembler")" || \
        fail "missing row for $sra:$srr:$read_type:$assembler in $summary_file"
    case "$note" in
        *"$expected"*) ;;
        *) fail "expected note for $sra:$srr:$read_type:$assembler to contain '$expected', saw: $note" ;;
    esac
}

write_summary_fixture() {
    # Write a tiny Nextflow fixture around the shared SUMMARY subworkflow.
    local fixture_file=$1

    cat > "$fixture_file" <<EOF
nextflow.enable.dsl=2

include { SUMMARY } from '${REPO_ROOT}/subworkflows/local/summary'

workflow {
  def outdir = file(params.outdir).toAbsolutePath().toString()
  def sample_summary = tuple(
    'main_sample',
    'main_sample',
    'ILLUMINA',
    'NovaSeq',
    'WGS',
    'short',
    'metaspades',
    file(params.summary_csv, checkIfExists: true)
  )

  def taxa_summary_ch = channel.of(sample_summary)
  def binning_note_entries_ch = channel.empty()

  if (params.fixture == 'multiple_notes') {
    binning_note_entries_ch = channel.of(
      tuple('main_sample', 'main_sample', 'ILLUMINA', 'NovaSeq', 'WGS', 'short', 'metaspades', 'metabat', file(params.metabat_note, checkIfExists: true)),
      tuple('main_sample', 'main_sample', 'ILLUMINA', 'NovaSeq', 'WGS', 'short', 'metaspades', 'dastool', file(params.dastool_note, checkIfExists: true))
    )
  }
  else if (params.fixture == 'empty_notes') {
    binning_note_entries_ch = channel.of(
      tuple('main_sample', 'main_sample', 'ILLUMINA', 'NovaSeq', 'WGS', 'short', 'metaspades', 'metabat', file(params.metabat_note, checkIfExists: true)),
      tuple('main_sample', 'main_sample', 'ILLUMINA', 'NovaSeq', 'WGS', 'short', 'metaspades', 'dastool', file(params.dastool_note, checkIfExists: true))
    )
  }

  SUMMARY(
    channel.empty(),
    channel.empty(),
    channel.empty(),
    channel.empty(),
    channel.empty(),
    channel.empty(),
    channel.empty(),
    channel.empty(),
    channel.empty(),
    taxa_summary_ch,
    binning_note_entries_ch,
    outdir
  )
}
EOF
}

run_nextflow() {
    # Run Nextflow with the repository helper scripts available on PATH.
    local log_file=$1
    shift

    if ! (cd "$TMP_ROOT" && PATH="$REPO_ROOT/bin:$PATH" nextflow run "$@" -ansi-log false > "$log_file" 2>&1); then
        sed -n '1,160p' "$log_file" >&2
        fail "Nextflow command failed; see $log_file"
    fi
}

run_summary_fixture() {
    # Run one SUMMARY fixture scenario and assert one main workflow row.
    local fixture_name=$1
    local expected_note_mode=$2
    local fixture_file="$TMP_ROOT/summary_fixture.nf"
    local outdir="$TMP_ROOT/summary_${fixture_name}_out"
    local workdir="$TMP_ROOT/summary_${fixture_name}_work"
    local log_file="$TMP_ROOT/summary_${fixture_name}.log"
    local summary_file="$outdir/summary.tsv"
    local metabat_note="$TMP_ROOT/${fixture_name}_metabat.note"
    local dastool_note="$TMP_ROOT/${fixture_name}_dastool.note"

    : > "$metabat_note"
    : > "$dastool_note"
    if [ "$fixture_name" = "multiple_notes" ]; then
        printf 'MetaBAT warning\n' > "$metabat_note"
        printf 'DAS Tool warning\n' > "$dastool_note"
    fi

    run_nextflow "$log_file" "$fixture_file" \
        --fixture "$fixture_name" \
        --summary_csv "$SUMMARY_CSV" \
        --metabat_note "$metabat_note" \
        --dastool_note "$dastool_note" \
        --binning \
        --outdir "$outdir" \
        -work-dir "$workdir"

    assert_file_exists "$summary_file"
    assert_summary_header "$summary_file"
    assert_data_row_count "$summary_file" 1
    assert_key_row_count "$summary_file" main_sample main_sample short metaspades 1

    if [ "$expected_note_mode" = "blank" ]; then
        assert_key_note_blank "$summary_file" main_sample main_sample short metaspades
    else
        assert_key_note_contains "$summary_file" main_sample main_sample short metaspades "metabat: MetaBAT warning"
        assert_key_note_contains "$summary_file" main_sample main_sample short metaspades "dastool: DAS Tool warning"
    fi
}

run_standalone_stub() {
    # Run standalone binning in stub mode and assert each sample appears once.
    local binning_tsv="$TMP_ROOT/binning_abs.tsv"
    local outdir="$TMP_ROOT/standalone_out"
    local workdir="$TMP_ROOT/standalone_work"
    local log_file="$TMP_ROOT/standalone.log"
    local summary_file="$outdir/summary.tsv"

    {
        printf 'sample\tread_type\treads\tassembly_fasta\n'
        printf 'short_sample\tshort\t%s,%s\t%s\n' \
            "$REPO_ROOT/test/binning/short_R1.fastq" \
            "$REPO_ROOT/test/binning/short_R2.fastq" \
            "$REPO_ROOT/test/binning/assembly.fasta"
        printf 'hifi_sample\thifi\t%s\t%s\n' \
            "$REPO_ROOT/test/binning/hifi.fastq" \
            "$REPO_ROOT/test/binning/assembly.fasta"
    } > "$binning_tsv"

    run_nextflow "$log_file" "$REPO_ROOT/binning.nf" \
        -stub-run \
        --binning_tsv "$binning_tsv" \
        --uniprot_db "$REPO_ROOT/test/binning/assembly.fasta" \
        --binners metabat \
        --refiners dastool \
        --outdir "$outdir" \
        -work-dir "$workdir"

    assert_file_exists "$summary_file"
    assert_summary_header "$summary_file"
    assert_data_row_count "$summary_file" 2
    assert_key_row_count "$summary_file" short_sample short_sample short provided 1
    assert_key_row_count "$summary_file" hifi_sample hifi_sample hifi provided 1
    assert_key_note_blank "$summary_file" short_sample short_sample short provided
    assert_key_note_blank "$summary_file" hifi_sample hifi_sample hifi provided
}

main() {
    # Run all binning summary row regression checks.
    TMP_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/binning_summary_rows_test.XXXXXX")"
    trap 'rm -rf "$TMP_ROOT"' EXIT

    SUMMARY_CSV="$TMP_ROOT/summary.csv"
    cat > "$SUMMARY_CSV" <<'EOF'
rank,ncbi_taxa,n_contigs,output_ids_csv,output_fasta
species,12345,7,ids.csv,records.fasta
EOF

    write_summary_fixture "$TMP_ROOT/summary_fixture.nf"

    run_summary_fixture no_notes blank
    run_summary_fixture multiple_notes noted
    run_summary_fixture empty_notes blank
    run_standalone_stub

    printf 'binning summary row tests passed\n'
}

main "$@"
