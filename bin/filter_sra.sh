#!/usr/bin/env bash
# Usage: ./filter_sra.sh SRR123456.metadata.tsv

in="$1"
out="${in%.metadata.tsv}.filtered.csv"

# Filters (edit easily)
PLATFORMS='^(ILLUMINA|BGISEQ|DNBSEQ|PACBIO_SMRT|OXFORD_NANOPORE)$'
SOURCES='^(GENOMIC|METAGENOMIC)$'
STRATEGIES='^WGS$'

awk -F'\t' -v OFS=',' -v P="$PLATFORMS" -v S="$SOURCES" -v T="$STRATEGIES" '
BEGIN {
  r="run_accession"; p="instrument_platform"; s="library_source"; t="library_strategy"
}
NR==1 {
  for(i=1;i<=NF;i++){ sub(/\r$/,"",$i); idx[$i]=i }
  print r, p, s, t; next
}
{
  sub(/\r$/,"",$0)
  if ($idx[p] ~ P && $idx[s] ~ S && $idx[t] ~ T)
    print $idx[r], $idx[p], $idx[s], $idx[t]
}
' "$in" > "$out"
