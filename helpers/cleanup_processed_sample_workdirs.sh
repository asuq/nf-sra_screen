#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd -P)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." >/dev/null 2>&1 && pwd -P)"
HELPER="$REPO_ROOT/external/nf-helper/helpers/cleanup_processed_sample_workdirs.sh"

exec "$HELPER" "$@"
