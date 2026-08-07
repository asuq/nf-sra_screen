#!/usr/bin/env bash

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${repo_root}"

export NXF_SYNTAX_PARSER=v2

profiles=(local slurm debug test oist gwdg marmic)
for profile in "${profiles[@]}"; do
  nextflow config -profile "${profile}" >/dev/null
done

echo "Validated Nextflow configuration for: ${profiles[*]}"
