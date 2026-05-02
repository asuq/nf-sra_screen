#!/usr/bin/env nextflow

/*
 * Return the command-line help text for the standalone binning workflow.
 */
def helpMessage() {
  log.info """
  ===============================
  nf-sra_screen: binning
  Standalone assembly-guided binning workflow
  Version: ${params.version}
  Author : Akito Shima (ASUQ)
  Email: asuq.4096@gmail.com
  ===============================
  Usage: nextflow run binning.nf [parameters]

  Required parameters:
    --binning_tsv   Path to binning.tsv
    --uniprot_db    Path to UniProt database (.dmnd)

  Optional parameters:
    --outdir        Output directory (default: ./output)
    --binners       Comma-separated binners (default: all-compatible; allowed: all-compatible,auto,metabat,semibin,rosella,comebin,vamb,lorbin)
    --refiners      Comma-separated refiners (default: binette; allowed: dastool,binette)
    --checkm2_db    CheckM2 DIAMOND database required with --refiners binette
    --semibin_environment  SemiBin2 pretrained environment (default: global)
    --gpu           Use GPU variants for COMEBin, VAMB, and HiFi-only LorBin
    --gpu_type      Optional GPU type for typed scheduler requests on GWDG
    --gpus          GPU count for scheduler requests on GWDG (default: 1)
    --max_retries   Maximum number of retries for each process (default: 3)
    --help          Show this help message

  binning.tsv columns:
    sample          Logical sample identifier
    read_type       short | nanopore | pacbio | hifi (local-read rows only)
    reads           Comma-separated FASTQ paths
    srr             Optional SRR accession to download raw reads
    assembly_fasta  Path to assembly FASTA

  Summary rows use read_type for the read class and assembler=provided.
  """.stripIndent()
}


include { STANDALONE_BINNING } from './subworkflows/local/standalone_binning'


/*
 * Abort with a consistent parameter error and show the help text.
 */
def missingParametersError() {
  log.error "Missing input parameters"
  helpMessage()
  error """
  For standalone binning, please provide:
    --binning_tsv and --uniprot_db
  """.stripIndent()
}



workflow {
    main:
    if (params.help) {
      helpMessage()
      exit 0
    }

    if (!params.binning_tsv || !params.uniprot_db) {
      missingParametersError()
    }

    STANDALONE_BINNING()

    workflow.onComplete = {
      def outdirPath = file(params.outdir ?: './output').toAbsolutePath()
      def summaryFile = outdirPath.resolve('summary.tsv')
      def traceFile = file("${workflow.launchDir}/execution-reports/trace.tsv").toAbsolutePath()
      def scriptFile = file("${workflow.projectDir}/bin/annotate_summary_from_trace.py").toAbsolutePath()

      log.info "onComplete: summary.tsv -> ${summaryFile}"
      log.info "onComplete: trace.tsv   -> ${traceFile}"
      log.info "onComplete: annotator   -> ${scriptFile}"

      if (!summaryFile.exists()) {
        log.warn "onComplete: ${summaryFile} not found; skipping scheduler annotation"
      }
      else if (!traceFile.exists()) {
        log.warn "onComplete: ${traceFile} not found; skipping scheduler annotation"
      }
      else {
        def cmd = [
          'python3',
          scriptFile.toString(),
          summaryFile.toString(),
          traceFile.toString()
        ]

        log.info "onComplete: running ${cmd.join(' ')}"

        def proc = new java.lang.ProcessBuilder(cmd)
          .directory(workflow.launchDir.toFile())
          .redirectError(java.lang.ProcessBuilder.Redirect.INHERIT)
          .redirectOutput(java.lang.ProcessBuilder.Redirect.INHERIT)
          .start()

        def rc = proc.waitFor()
        if (rc != 0) {
          log.warn "onComplete: annotator script exited with code ${rc}"
        }
        else {
          log.info "onComplete: summary.tsv successfully annotated with scheduler error information"
        }
      }
    }
}
