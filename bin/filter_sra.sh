#!/usr/bin/env bash
# Usage: ./filter_sra.sh ACCESSION.metadata.tsv
# Output: ACCESSION.filtered.csv
# Columns: accession,run_accession,instrument_platform,instrument_model,library_source,library_strategy,assembler

in="$1"
ACCESSION="${in%.metadata.tsv}"
out="${ACCESSION}.filtered.csv"

# Filters
PLATFORMS='^(ILLUMINA|BGISEQ|DNBSEQ|PACBIO_SMRT|OXFORD_NANOPORE)$'
SOURCES='^(GENOMIC|METAGENOMIC)$'
STRATEGIES='^WGS$'

awk -F'\t' -v OFS=',' \
    -v P="$PLATFORMS" -v S="$SOURCES" -v T="$STRATEGIES" \
    -v A="$ACCESSION" '
BEGIN {
  r="run_accession"
  p="instrument_platform"
  m="instrument_model"
  s="library_source"
  t="library_strategy"
  print "accession", r, p, m, s, t, "assembler"
}
NR==1 {
  for (i=1;i<=NF;i++) { sub(/\r$/,"",$i); idx[$i]=i }
  next
}
{
  sub(/\r$/,"",$0)
  if (!(r in idx) || !(p in idx) || !(m in idx) || !(s in idx) || !(t in idx)) next

  plat  = $idx[p]
  src   = $idx[s]
  strat = $idx[t]

  # Apply whitelist filters
  if (plat ~ P && src ~ S && strat ~ T) {
    run   = $idx[r]
    model = $idx[m]

    # normalize model for exact matching:
    lmodel = tolower(model)
    gsub(/[[:space:]]+/, " ", lmodel)                    # collapse spaces
    sub(/^[[:space:]]+/, "", lmodel); sub(/[[:space:]]+$/, "", lmodel)  # trim

    # Decide assembler:
    # - ILLUMINA/DNBSEQ/BGISEQ -> short
    # - OXFORD_NANOPORE -> long_nano
    # - PACBIO_SMRT with Pacbio RS/RSII, Sequel in model -> long_pacbio
    # - PACBIO_SMRT other models (e.g., Sequel II, Sequel IIe, Revio) -> long_hifi

    asm = ""
    if (plat ~ /^(ILLUMINA|DNBSEQ|BGISEQ)$/) {
      asm = "short"
    } else if (plat == "OXFORD_NANOPORE") {
      asm = "long_nano"
    } else if (plat == "PACBIO_SMRT") {
      # exact models for non-HiFi long reads:
      if (lmodel == "pacbio rs" || lmodel == "pacbio rs ii" || lmodel == "sequel") {
        asm = "long_pacbio"
      } else {
        # everything else under PACBIO_SMRT → HiFi bucket (e.g., "sequel ii", "sequel iie", "revio", etc.)
        asm = "long_hifi"
      }
    }

    print A, run, plat, model, src, strat, asm
  }
}
' "$in" > "$out"
