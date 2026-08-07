#!/usr/bin/env bash

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." >/dev/null 2>&1 && pwd -P)"

fail() {
  printf 'FAIL: %s\n' "$*" >&2
  exit 1
}

run_selection() {
  NXF_OFFLINE=true nextflow -C /dev/null run "$REPO_ROOT/tests/assembler_selection.nf" \
    -ansi-log false \
    "$@" \
    2>&1
}

assert_selection() {
  local output=$1
  local read_type=$2
  local expected=$3
  local line

  printf -v line 'ASSEMBLERS\t%s\t%s' "$read_type" "$expected"
  grep -Fqx -- "$line" <<<"$output" \
    || fail "expected selection '$line' in Nextflow output"
}

auto_output="$(run_selection --assembler auto)"
assert_selection "$auto_output" short metaspades
assert_selection "$auto_output" nanopore metaflye
assert_selection "$auto_output" pacbio metaflye
assert_selection "$auto_output" hifi myloasm

myloasm_output="$(run_selection --assembler myloasm)"
assert_selection "$myloasm_output" short ""
assert_selection "$myloasm_output" nanopore myloasm
assert_selection "$myloasm_output" pacbio ""
assert_selection "$myloasm_output" hifi myloasm

plural_myloasm_output="$(run_selection --assemblers myloasm)"
assert_selection "$plural_myloasm_output" nanopore myloasm
assert_selection "$plural_myloasm_output" hifi myloasm

metaflye_output="$(run_selection --assembler metaflye)"
assert_selection "$metaflye_output" short ""
assert_selection "$metaflye_output" nanopore metaflye
assert_selection "$metaflye_output" pacbio metaflye
assert_selection "$metaflye_output" hifi metaflye

all_output="$(run_selection --assembler all)"
assert_selection "$all_output" short metaspades,unicycler
assert_selection "$all_output" nanopore metaflye,myloasm
assert_selection "$all_output" pacbio metaflye
assert_selection "$all_output" hifi metaflye,myloasm

flye_alias_output="$(run_selection --assembler flye)"
assert_selection "$flye_alias_output" nanopore metaflye

printf 'Assembler selection tests passed.\n'
