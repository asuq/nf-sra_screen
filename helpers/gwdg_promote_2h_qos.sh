#!/usr/bin/env bash
#
# Promote eligible pending GWDG jobs into the 2h QOS when free slots exist.
#
# This helper is intentionally separate from Nextflow config. Run it from a
# login-node tmux session while the pipeline is active if you want short pending
# jobs to opportunistically use the small 2h QOS capacity.

set -euo pipefail

CAP=10
INTERVAL_SECONDS=60
QUIET=0
ONCE=0

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Options:
  --cap N              Maximum number of jobs allowed in QOS=2h (default: 10)
  --interval SECONDS   Seconds between checks (default: 60)
  --once               Run one check and exit
  -q, --quiet          Hide routine status lines, but still print job updates
  -h, --help           Show this help message

The helper promotes pending jobs owned by USER when:
  - their current QOS is not 2h
  - their requested walltime is at most 2 hours
  - free 2h QOS slots are available
EOF
}

fail() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

require_command() {
    command -v "$1" >/dev/null 2>&1 || fail "required command not found: $1"
}

require_positive_integer() {
    local name=$1
    local value=$2

    case "$value" in
        ''|*[!0-9]*)
            fail "$name must be a positive integer, got '$value'"
            ;;
    esac

    if [ "$value" -le 0 ]; then
        fail "$name must be greater than zero, got '$value'"
    fi
}

log_status() {
    if [ "$QUIET" -eq 0 ]; then
        printf '%s  current_2h=%s  free_slots=%s\n' \
            "$(date '+%F %T')" "$1" "$2"
    fi
}

promote_once() {
    local current
    local free

    current=$(squeue -h -u "$USER" -q 2h -o "%i" | wc -l | awk '{ print $1 }')
    free=$((CAP - current))

    log_status "$current" "$free"

    if [ "$free" -le 0 ]; then
        return 0
    fi

    squeue -h -u "$USER" -t PD -o "%i|%q|%l|%j" \
        | awk -F'|' -v limit="$free" '
            function tsec(t, a, b, n, d) {
                d = 0

                if (t == "UNLIMITED" || t == "NOT_SET" || t == "N/A") {
                    return 999999999
                }

                if (index(t, "-")) {
                    split(t, b, "-")
                    d = b[1]
                    t = b[2]
                }

                n = split(t, a, ":")

                if (n == 3) return d * 86400 + a[1] * 3600 + a[2] * 60 + a[3]
                if (n == 2) return d * 86400 + a[1] * 60 + a[2]
                if (n == 1) return d * 86400 + a[1]

                return 999999999
            }

            $2 != "2h" && tsec($3) <= 7200 {
                print $1
                promoted++
                if (promoted >= limit) {
                    exit
                }
            }
        ' \
        | while IFS= read -r job_id; do
            if [ -n "$job_id" ]; then
                printf 'Changing JobId=%s to QOS=2h\n' "$job_id"
                scontrol update JobId="$job_id" QOS=2h
            fi
        done
}

while [ "$#" -gt 0 ]; do
    case "$1" in
        --cap)
            [ "$#" -ge 2 ] || fail "--cap requires a value"
            CAP=$2
            shift 2
            ;;
        --cap=*)
            CAP=${1#*=}
            shift
            ;;
        --interval)
            [ "$#" -ge 2 ] || fail "--interval requires a value"
            INTERVAL_SECONDS=$2
            shift 2
            ;;
        --interval=*)
            INTERVAL_SECONDS=${1#*=}
            shift
            ;;
        --once)
            ONCE=1
            shift
            ;;
        -q|--quiet)
            QUIET=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            fail "unknown option: $1"
            ;;
    esac
done

: "${USER:?USER must be set}"

require_positive_integer "--cap" "$CAP"
require_positive_integer "--interval" "$INTERVAL_SECONDS"
require_command awk
require_command date
require_command scontrol
require_command squeue
require_command wc

while :; do
    promote_once

    if [ "$ONCE" -eq 1 ]; then
        break
    fi

    sleep "$INTERVAL_SECONDS"
done
