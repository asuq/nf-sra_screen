#!/usr/bin/env bash
# Archive one raw binner directory and remove the unpacked copy.
#
# Args:
#   --dir      Raw binner directory to archive.
#   --archive  Output archive path. Defaults to <dir>.tar.gz.

set -euo pipefail

bin_dir=""
archive=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dir) bin_dir="$2"; shift 2 ;;
    --archive) archive="$2"; shift 2 ;;
    *)
      echo "archive_binner_dir.sh: unknown argument: $1" >&2
      exit 1 ;;
  esac
done

if [[ -z "$bin_dir" ]]; then
  echo "archive_binner_dir.sh: --dir is required" >&2
  exit 1
fi

if [[ ! -d "$bin_dir" ]]; then
  echo "archive_binner_dir.sh: directory does not exist: ${bin_dir}" >&2
  exit 1
fi

if [[ -z "$archive" ]]; then
  archive="${bin_dir}.tar.gz"
fi

bin_parent="$(dirname "$bin_dir")"
bin_name="$(basename "$bin_dir")"

if [[ -z "$bin_name" || "$bin_name" == "." || "$bin_name" == "/" ]]; then
  echo "archive_binner_dir.sh: refusing unsafe directory: ${bin_dir}" >&2
  exit 1
fi

archive_tmp="${archive}.tmp.$$"

rm -f "$archive" "$archive_tmp"
tar -czf "$archive_tmp" -C "$bin_parent" "$bin_name"
mv "$archive_tmp" "$archive"
rm -rf "$bin_dir"
