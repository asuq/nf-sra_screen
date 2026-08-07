#!/usr/bin/env bash

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${repo_root}"

export NXF_SYNTAX_PARSER=v2

fail() {
  printf 'FAIL: %s\n' "$*" >&2
  exit 1
}

profiles=(local slurm debug test oist gwdg marmic)
for profile in "${profiles[@]}"; do
  nextflow config -profile "${profile}" >/dev/null
done

oist_config="$(nextflow config -profile oist -flat)"
grep -F -- "singularity.enabled = false" <<<"${oist_config}" >/dev/null \
  || fail "OIST profile must disable Singularity"
grep -F -- "apptainer.enabled = true" <<<"${oist_config}" >/dev/null \
  || fail "OIST profile must enable Apptainer"

nextflow inspect -profile oist -ignore-errors . \
  | python -c '
import json
import sys

expected = "quay.io/asuq1617/iseq:1.9.8-sratools3.4.1-r1@sha256:59e96013353dbecddf483aac16591c591b429219b7857a9b277c3bfcb1060ed1"
inspection = json.load(sys.stdin)
downloads = [
    process
    for process in inspection["processes"]
    if process.get("name") == "DOWNLOAD_SRR"
]
if len(downloads) != 1:
    raise SystemExit(f"expected one DOWNLOAD_SRR process, found {len(downloads)}")
actual = downloads[0].get("container")
if actual != expected:
    raise SystemExit(f"DOWNLOAD_SRR container mismatch: expected {expected!r}, found {actual!r}")
'

echo "Validated Nextflow configuration for: ${profiles[*]}"
