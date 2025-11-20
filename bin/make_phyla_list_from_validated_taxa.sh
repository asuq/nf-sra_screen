#!/usr/bin/env bash
# Usage: make_phyla_list_from_validated_taxa.sh validated_taxa.csv phyla_to_check.txt

set -euo pipefail

in_csv="$1"
out_txt="$2"

awk -F',' '
  NR==1 {
    for (i = 1; i <= NF; i++) if ($i == "gtdb_phylum") { c = i; break }
    next
  }
  c && $c != "" { seen[$c] = 1 }

  END {
    print "phyla"
    for (v in seen) print v
  }
' "$in_csv" > "$out_txt"
