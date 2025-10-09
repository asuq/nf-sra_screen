#!/usr/bin/env nextflow
nextflow.enable.dsl=2


//-- Help Message ---------------------------------------------------------------

def helpMessage() {
  log.info """
  ===============================
  nf-sra_screen
  Nextflow pipeline for screening SRA genomes
  Version: 0.1.0
  Author : Akito Shima (ASUQ)
  Email: asuq.4096@gmail.com
  ===============================
  Usage: nextflow run main.nf [parameters]

  Required parameters:
    --sra           Path to sra.csv (header: sra)
    --nt_db         Path to BLAST nt database
    --blobtools_db  Path to BlobTools database

  Optional parameters:
    --help          Show this help message
    --outdir        Output directory (default: ./output)
    --spades_opts   SPAdes options (default: --threads \$task.cpus)
  """.stripIndent()
}


def missingParametersError() {
    log.error "Missing input parameters"
    helpMessage()
    error "Please provide all required parameters: --sra, --nt_db, and --blobtools_db"
}


//-- Processes -----------------------------------------------------------------

// Download SRA
process download_sra_metadata {
    tag "${sra}"
    publishDir "${params.outdir}/metadata/${sra}", mode: 'copy', overwrite: true

    input:
    val sra

    output:
    tuple val(sra), path("${sra}.filtered.csv"), emit: filtered_csv

    script:
    """
    # Retrieve metadata
    iseq -m -i "${sra}"

    # Extract SRR to screen
    filter_sra.sh "${sra}.metadata.tsv"
    """
}


//-- Workflow ------------------------------------------------------------------
workflow {
  // Parameter parsing
  if (params.help) {
    helpMessage()
    exit 0
  }

  //TODO: Add check for nt_db and blobtools_db
  if (!params.sra) {
    missingParametersError()
    exit 1
  }

  // Channel setup
  sra_ch = Channel.fromPath(params.sra, checkIfExists: true)
                  .splitCsv(header: true, strip: true)
                  .map { row -> row.sra.trim() }
                  .filter { it }
                  .distinct()

  // Step 1: extract metadata & filter SRR
  filtered_srr = download_sra_metadata(sra_ch).view()
}
