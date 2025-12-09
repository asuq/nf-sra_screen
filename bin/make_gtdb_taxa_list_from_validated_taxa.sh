#!/usr/bin/env bash
# Usage: make_gtdb_taxa_list_from_validated_taxa.sh validated_taxa.csv taxa_to_check.txt
#
# Semantics:
#   - validated_taxa.csv comes from validate_taxa.py and has columns:
#       rank,taxa,ncbi_phylum,gtdb_phylum
#   - For Sandpiper / SingleM filtering:
#       * If 'taxa' is GTDB style (d__/p__/c__/o__/f__/g__/s__), use that value.
#       * Otherwise, use the value from 'gtdb_phylum' (if non-empty).
#   - Output is a text file with header 'gtdb_taxa', one target per line.

set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 validated_taxa.csv taxa_to_check.txt" >&2
  exit 1
fi

in_csv="$1"
out_txt="$2"

if [[ ! -f "$in_csv" ]]; then
  echo "ERROR: input CSV not found: ${in_csv}" >&2
  exit 1
fi


# header: rank,taxa,ncbi_phylum,gtdb_phylum
{
  echo "gtdb_taxa"

  awk -F',' '
    BEGIN {
      col_taxa = 0
      col_gtdb = 0
    }

    NR == 1 {
      # Header line: detect column indices (case-insensitive).
      for (i = 1; i <= NF; i++) {
        h = $i
        gsub(/^[ \t"]+|[ \t"]+$/, "", h)
        lh = tolower(h)
        if (lh == "taxa") {
          col_taxa = i
        } else if (lh == "gtdb_phylum" || lh == "gtdb_taxa") {
          col_gtdb = i
        }
      }
      if (!col_taxa) {
        print "ERROR: header missing a 'taxa' column in " FILENAME > "/dev/stderr"
        exit 1
      }
      # gtdb_phylum is optional; we just won’t use it if absent.
      next
    }

    NR > 1 {
      if (!col_taxa) {
        # Should not happen if header was parsed correctly.
        next
      }

      taxa = $col_taxa
      gtdb = (col_gtdb ? $col_gtdb : "")

      # Strip outer quotes and whitespace.
      gsub(/^"|"$/, "", taxa)
      gsub(/^"|"$/, "", gtdb)
      sub(/^[ \t]+/, "", taxa); sub(/[ \t]+$/, "", taxa)
      sub(/^[ \t]+/, "", gtdb); sub(/[ \t]+$/, "", gtdb)

      if (taxa == "" && gtdb == "") {
        next
      }

      # GTDB style detection (case-insensitive): d__/p__/c__/o__/f__/g__/s__
      ltaxa = tolower(taxa)
      is_gtdb = (ltaxa ~ /^(d__|p__|c__|o__|f__|g__|s__)/)

      target = ""
      if (is_gtdb) {
        # Use the GTDB-style taxon exactly as given (e.g. p__Bacillota).
        target = taxa
      } else {
        # NCBI input: use mapped GTDB phylum if available.
        if (gtdb == "" || tolower(gtdb) == "nan" || tolower(gtdb) == "none") {
          # No GTDB phylum mapping (e.g. Eukaryotes, Viruses, or unmapped).
          next
        }
        # User requested to "use the value of gtdb_phylum" directly.
        target = gtdb
      }

      if (target == "") {
        next
      }

      # Deduplicate by raw target string.
      if (!(target in seen)) {
        seen[target] = 1
        print target
      }
    }
  ' "$in_csv"
} > "$out_txt"
