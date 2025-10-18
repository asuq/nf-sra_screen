#!/bin/bash
# Usage: ./filter_sra.sh ACCESSION.metadata.(tsv|csv)
# Output:
#   ACCESSION.filtered.csv  # kept SRRs
#   ACCESSION.skipped.csv   # skipped SRRs with reason
#
# Columns (filtered):
#   accession,run_accession,instrument_platform,instrument_model,library_source,library_strategy,assembler
#
# Columns (skipped):
#   accession,run_accession,instrument_platform,instrument_model,library_source,library_strategy,skip_reason

set -euo pipefail

acc="${1:-}"
if [[ -z "$acc" ]]; then
  echo "ERROR: missing accession. Usage: $0 ACCESSION" >&2
  exit 1
fi

# Derive ACCESSION prefix for either .metadata.tsv or .metadata.csv
in=""
if   [[ -f "${acc}.metadata.tsv" ]]; then in="${acc}.metadata.tsv"
elif [[ -f "${acc}.metadata.csv" ]]; then in="${acc}.metadata.csv"
else
  echo "ERROR: metadata not found: expected '${acc}.metadata.tsv' or '${acc}.metadata.csv' in current directory." >&2
  exit 1
fi

out_kept="${acc}.filtered.csv"
out_skip="${acc}.skipped.csv"

# Whitelist filters (edit easily)
ACCESSION="$acc"
PLATFORMS='^(ILLUMINA|BGISEQ|DNBSEQ|PACBIO_SMRT|OXFORD_NANOPORE)$'
SOURCES='^(GENOMIC|METAGENOMIC)$'
STRATEGIES='^WGS$'

# Always create headers
printf "accession,run_accession,instrument_platform,instrument_model,library_source,library_strategy,assembler\n" > "$out_kept"
printf "accession,run_accession,instrument_platform,instrument_model,library_source,library_strategy,skip_reason\n" > "$out_skip"

header="$(head -n1 "$in")"

# ------------------------------ GSA CSV branch -------------------------------
# GSA header starts with: Run,Center,ReleaseDate,FileType,FileName,FileSize,Download_path,...
if [[ "$header" == Run,Center,ReleaseDate,FileType,FileName,FileSize,Download_path,* ]]; then
  awk -F',' -v OFS=',' \
      -v P="$PLATFORMS" -v S="$SOURCES" -v T="$STRATEGIES" \
      -v A="$ACCESSION" -v OUT_KEEP="$out_kept" -v OUT_SKIP="$out_skip" '
  function trim(x){ sub(/^[[:space:]]+/,"",x); sub(/[[:space:]]+$/,"",x); return x }
  NR==1{
    for(i=1;i<=NF;i++){ gsub(/\r$/,"",$i); idx[$i]=i }
    # Required GSA columns
    needR="Run"; needM="Platform"; needS="LibrarySource"; needT="LibraryStrategy"
    if(!(needR in idx) || !(needM in idx) || !(needS in idx) || !(needT in idx)){
      print "ERROR: Missing required GSA columns (Run, Platform, LibrarySource, LibraryStrategy)." > "/dev/stderr"
      exit 1
    }
    ir=idx[needR]; im=idx[needM]; is=idx[needS]; it=idx[needT]
    next
  }
  {
    gsub(/\r$/,"")
    run=$ir; model=$im; src=$is; strat=$it

    # Derive instrument_platform from the model string (GSA puts model in "Platform")
    lm=tolower(model); gsub(/[[:space:]]+/, " ", lm); lm=trim(lm)
    plat=""
    if     (lm ~ /illumina/)                                   plat="ILLUMINA"
    else if(lm ~ /dnbseq/)                                     plat="DNBSEQ"
    else if(lm ~ /mgiseq|bgiseq|(^|[^a-z])bgi($|[^a-z])|(^|[^a-z])mgi($|[^a-z])/) plat="BGISEQ"
    else if(lm ~ /oxford|nanopore|minion|gridion|promethion|\bont\b/)             plat="OXFORD_NANOPORE"
    else if(lm ~ /pacbio|sequel|revio|\brs\b/)                 plat="PACBIO_SMRT"

    keep = (plat ~ P) && (src ~ S) && (strat ~ T)

    if (keep){
      # Assembler decision (based on derived platform + model keywords)
      asm=""
      if (plat ~ /^(ILLUMINA|DNBSEQ|BGISEQ)$/) asm="short"
      else if (plat=="OXFORD_NANOPORE")        asm="long_nano"
      else if (plat=="PACBIO_SMRT"){
        if (lm ~ /(^|[^a-z])rs($|[^a-z])|(^|[^a-z])rs[[:space:]]*ii($|[^a-z])|(^|[^a-z])sequel($|[^a-z])/ &&
            lm !~ /sequel[[:space:]]*ii|sequel[[:space:]]*2|iie|revio/) asm="long_pacbio"
        else asm="long_hifi"
      } else asm="short"

      print A, run, plat, model, src, strat, asm >> OUT_KEEP
    } else {
      reason = (src !~ S ? "source" : (strat !~ T ? "strategy" : (plat !~ P ? "platform" : "other")))
      print A, run, plat, model, src, strat, reason >> OUT_SKIP
    }
  }' "$in"
  exit 0
fi

# --------------------- ENA TSV / NCBI-TSV generic branch ---------------------
# Tab-delimited cases (ENA fields=all or NCBI RunInfo TSV via iSeq)
awk -F'\t' -v OFS=',' \
    -v P="$PLATFORMS" -v S="$SOURCES" -v T="$STRATEGIES" \
    -v A="$ACCESSION" -v OUT_KEEP="$out_kept" -v OUT_SKIP="$out_skip" '
function trim(x){ sub(/^[[:space:]]+/,"",x); sub(/[[:space:]]+$/,"",x); return x }
function normkey(x,   y){ y=tolower(x); gsub(/[^a-z0-9]+/,"_",y); return y }
function pick(arr,  i) { for (i=1; i in arr; i++) if (arr[i] in idx) return arr[i]; return "" }

BEGIN{
  want_run="run_accession"; want_plat="instrument_platform";
  want_model="instrument_model"; want_src="library_source"; want_strat="library_strategy";
}
NR==1{
  for(i=1;i<=NF;i++){
    sub(/\r$/,"",$i)
    k = normkey($i)
    idx[k] = i
  }

  # Aliases to cover ENA + NCBI-TSV
  ralts[1]="run_accession"; ralts[2]="run"
  palts[1]="instrument_platform"; palts[2]="platform"
  malts[1]="instrument_model"; malts[2]="model"
  salts[1]="library_source"; salts[2]="librarysource"
  talts[1]="library_strategy"; talts[2]="librarystrategy"

  rk = pick(ralts); pk = pick(palts); mk = pick(malts); sk = pick(salts); tk = pick(talts)

  if (rk=="" || pk=="" || mk=="" || sk=="" || tk=="") {
    print "ERROR: Required columns missing after header normalization." > "/dev/stderr"
    print "Looked for any of:" > "/dev/stderr"
    print "  run: run_accession|run" > "/dev/stderr"
    print "  platform: instrument_platform|platform" > "/dev/stderr"
    print "  model: instrument_model|instrumentmodel|model" > "/dev/stderr"
    print "  source: library_source|librarysource" > "/dev/stderr"
    print "  strategy: library_strategy|librarystrategy" > "/dev/stderr"
    exit 1
  }

  ir = idx[rk]; ip = idx[pk]; im = idx[mk]; is = idx[sk]; it = idx[tk]
  next
}
{
  sub(/\r$/,"",$0)
  run   = $(ir)
  plat  = $(ip)
  model = $(im)
  src   = $(is)
  strat = $(it)

  keep = (plat ~ P) && (src ~ S) && (strat ~ T)

  if (keep) {
    lmodel = tolower(model); gsub(/[[:space:]]+/, " ", lmodel); lmodel = trim(lmodel)

    asm = ""
    if (plat ~ /^(ILLUMINA|DNBSEQ|BGISEQ)$/) {
      asm = "short"
    } else if (plat == "OXFORD_NANOPORE") {
      asm = "long_nano"
    } else if (plat == "PACBIO_SMRT") {
      if (lmodel ~ /(^|[^a-z])rs($|[^a-z])|(^|[^a-z])rs[[:space:]]*ii($|[^a-z])|(^|[^a-z])sequel($|[^a-z])/ &&
          lmodel !~ /sequel[[:space:]]*ii|sequel[[:space:]]*2|iie|revio/) asm = "long_pacbio"
      else asm = "long_hifi"
    } else {
      asm = "short"
    }

    print A, run, plat, model, src, strat, asm >> OUT_KEEP

  } else {
    reason = (src !~ S ? "source" : (strat !~ T ? "strategy" : (plat !~ P ? "platform" : "other")))
    print A, run, plat, model, src, strat, reason >> OUT_SKIP
  }
}
' "$in"
