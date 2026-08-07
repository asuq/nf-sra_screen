#!/usr/bin/env bash

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." >/dev/null 2>&1 && pwd -P)"
TMP_ROOT="$(mktemp -d /tmp/nf-sra-screen-myloasm-test.XXXXXX)"
trap 'rm -rf -- "$TMP_ROOT"' EXIT

fail() {
  printf 'FAIL: %s\n' "$*" >&2
  exit 1
}

assert_file_exists() {
  [[ -f "$1" ]] || fail "expected file to exist: $1"
}

assert_file_missing() {
  [[ ! -e "$1" ]] || fail "expected file to be absent: $1"
}

assert_file_contains() {
  local file=$1
  local expected=$2

  grep -Fq -- "$expected" "$file" || fail "expected '$expected' in $file"
}

write_fake_tools() {
  local fake_bin=$1

  mkdir -p "$fake_bin"

  cat > "$fake_bin/myloasm" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

printf '%s\n' "$@" > "${FAKE_MYLOASM_ARGS}"
printf '>u1ctg_len-4_circular-no_depth-1-1-1_duplicated-no\nACGT\n' > assembly_primary.fa
printf 'H\tVN:Z:1.0\nS\tu1ctg\t*\tLN:i:3\n' > final_contig_graph.gfa
EOF

  cat > "$fake_bin/mylotools" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

printf '%s\n' "$@" > "${FAKE_MYLOTOOLS_ARGS}"
[[ "$1" == "annotate-gfa" ]]
shift

gfa=""
fasta=""
output=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --gfa) gfa=$2; shift 2 ;;
    --fasta) fasta=$2; shift 2 ;;
    --output) output=$2; shift 2 ;;
    *) exit 64 ;;
  esac
done

[[ "$gfa" == "final_contig_graph.gfa" ]]
[[ "$fasta" == "assembly.fasta" ]]
[[ "$output" == "assembly.gfa" ]]
[[ -f "$gfa" && -f "$fasta" ]]
[[ ! -e assembly_primary.fa ]]

printf 'partial graph\n' > "$output"
if [[ "${FAKE_MYLOTOOLS_EXIT:-0}" != "0" ]]; then
  exit "${FAKE_MYLOTOOLS_EXIT}"
fi
printf 'annotated from %s using %s\n' "$fasta" "$gfa" > "$output"
EOF

  chmod +x "$fake_bin/myloasm" "$fake_bin/mylotools"
}

run_wrapper() {
  local work_dir=$1
  local fake_bin=$2
  local attempt=$3
  local max_retries=$4
  local mylotools_exit=$5

  mkdir -p "$work_dir"
  printf '@read\nACGT\n+\nIIII\n' > "$work_dir/reads.fastq"

  (
    cd "$work_dir"
    FAKE_MYLOASM_ARGS="$work_dir/myloasm.args" \
    FAKE_MYLOTOOLS_ARGS="$work_dir/mylotools.args" \
    FAKE_MYLOTOOLS_EXIT="$mylotools_exit" \
    PATH="$fake_bin:$REPO_ROOT/bin:$PATH" \
      run_myloasm_hifi.sh \
        --reads reads.fastq \
        --cpus 3 \
        --attempt "$attempt" \
        --max-retries "$max_retries"
  )
}

test_successful_annotation() {
  local fake_bin="$TMP_ROOT/fake-success-bin"
  local work_dir="$TMP_ROOT/success"

  write_fake_tools "$fake_bin"
  run_wrapper "$work_dir" "$fake_bin" 1 2 0

  assert_file_exists "$work_dir/assembly.fasta"
  assert_file_exists "$work_dir/assembly.gfa"
  assert_file_contains "$work_dir/assembly.gfa" \
    "annotated from assembly.fasta using final_contig_graph.gfa"
  assert_file_missing "$work_dir/assembly_primary.fa"
  assert_file_missing "$work_dir/FAIL.note"

  diff -u <(printf '%s\n' \
    annotate-gfa \
    --gfa final_contig_graph.gfa \
    --fasta assembly.fasta \
    --output assembly.gfa) \
    "$work_dir/mylotools.args"
}

test_annotation_retry_failure() {
  local fake_bin="$TMP_ROOT/fake-retry-bin"
  local work_dir="$TMP_ROOT/retry"

  write_fake_tools "$fake_bin"
  if run_wrapper "$work_dir" "$fake_bin" 1 2 42 > "$work_dir.log" 2>&1; then
    fail "annotation retry attempt unexpectedly succeeded"
  fi

  assert_file_exists "$work_dir/assembly.fasta"
  assert_file_missing "$work_dir/assembly.gfa"
  assert_file_missing "$work_dir/FAIL.note"
  assert_file_contains "$work_dir.log" "myloasm: GFA annotation failed"
}

test_annotation_final_soft_failure() {
  local fake_bin="$TMP_ROOT/fake-final-bin"
  local work_dir="$TMP_ROOT/final"

  write_fake_tools "$fake_bin"
  if ! run_wrapper "$work_dir" "$fake_bin" 2 2 42 > "$work_dir.log" 2>&1; then
    sed -n '1,120p' "$work_dir.log" >&2
    fail "final annotation soft-failure attempt failed"
  fi

  assert_file_exists "$work_dir/assembly.fasta"
  assert_file_missing "$work_dir/assembly.gfa"
  assert_file_contains "$work_dir/FAIL.note" "myloasm: GFA annotation failed"
}

test_successful_annotation
test_annotation_retry_failure
test_annotation_final_soft_failure

printf 'Myloasm GFA annotation wrapper tests passed.\n'
