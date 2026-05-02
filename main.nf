#!/usr/bin/env nextflow

//-- Help Message ---------------------------------------------------------------

def helpMessage() {
  log.info """
  ===============================
  nf-sra_screen
  Nextflow pipeline for screening SRA genomes
  Version: ${params.version}
  Author : Akito Shima (ASUQ)
  Email: asuq.4096@gmail.com
  ===============================
  Usage: nextflow run main.nf [parameters]

  Required parameters (assembly + binning):
    (A) SRA mode:
    --sra           Path to sra.csv (header: sra)
    (B) FASTQ mode:
    --fastq_tsv     Path to fastq.tsv (header: sample,read_type,reads)
                    (reads: comma-separated FASTQ paths)
    You can combine both modes by providing both --sra and --fastq_tsv

    And:
    --taxdump       Path to taxdump database folder
    --uniprot_db    Path to Uniprot database (.dmnd)

  Optional parameters (enable target taxa screening & extraction):
    --taxa          Path to taxa.csv for extraction (header: rank, taxa)
    --gtdb_ncbi_map Path to folder with GTDB-NCBI mapping Excel files
    --sandpiper_db  Path to Sandpiper database folder
    --singlem_db    Path to SingleM database folder

  Misc:
    --help          Show this help message
    --noassembly    Skip ASSEMBLY/BINNING
    --binning       Also run BINNING after ASSEMBLY
    --assemblers    Assembly tools: auto, all, or comma-separated tool names
                    (default: auto; aliases: --assembler, spades, flye)
    --binners       Comma-separated binners (default: all-compatible; allowed: all-compatible,auto,metabat,semibin,rosella,comebin,vamb,lorbin)
    --refiners      Comma-separated refiners (default: binette; allowed: dastool,binette)
    --checkm2_db    CheckM2 DIAMOND database required with --refiners binette
    --semibin_environment  SemiBin2 pretrained environment (default: global)
    --gpu           Use GPU variants for COMEBin, VAMB, and HiFi-only LorBin
    --gpu_type      Optional GPU type for typed scheduler requests on GWDG
    --gpus          GPU count for scheduler requests on GWDG (default: 1)
    --outdir        Output directory (default: ./output)
    --max_retries   Maximum number of retries for each process (default: 3)
  """.stripIndent()
}


def missingParametersError() {
  log.error "Missing input parameters"
  helpMessage()

  def noAssembly = params.noassembly?.toString()?.toBoolean() ?: false

  if (noAssembly) {
    error """
    For --noassembly (only PRE_SCREENING), please provide:
      at least one of --sra or --fastq_tsv
    If you enable taxa screening (--taxa), please also provide:
      --taxdump, --gtdb_ncbi_map, and --singlem_db
      and for sra mode also --sandpiper_db
    """.stripIndent()
  }
  else {
    error """
    For assembly + binning, please provide:
      --taxdump and --uniprot_db
      and at least one of --sra or --fastq_tsv
    If you also want to enable target taxa screening & extraction, please provide:
      --taxa, --gtdb_ncbi_map, and --singlem_db
      and for sra mode also --sandpiper_db
    """.stripIndent()
  }
}


include {
  selectedAssemblerTokens
} from './lib/workflow_helpers.nf'


include { VALIDATE_TAXA } from './modules/local/validate_taxa'
include { PRE_SCREENING } from './subworkflows/local/pre_screening'
include { FASTQ_PRE_SCREENING } from './subworkflows/local/fastq_pre_screening'

include { ASSEMBLY } from './subworkflows/local/assembly'

include { BINNING } from './subworkflows/local/binning'

include { CREATE_EMPTY_SUMMARY } from './modules/local/create_empty_summary'
include { SUMMARY } from './subworkflows/local/summary'

workflow {
    main:
    if (params.help) {
      helpMessage()
      exit 0
    }

    if (params.gpu?.toString()?.toBoolean()) {
      error "GPU mode is planned but not implemented yet for the current binning selection workflow"
    }

    def sraMode   = params.sra != null
    def fastqMode = params.fastq_tsv != null

    if (!sraMode && !fastqMode) {
      log.error "Error: either --sra or --fastq_tsv parameter must be provided"
      missingParametersError()
    }

    def doScreening = params.taxa != null
    def noAssembly = params.noassembly?.toString()?.toBoolean() ?: false
    def doAssembly = !noAssembly
    if (doAssembly) {
      selectedAssemblerTokens()
    }

    // If noassembly, binning makes no sense.
    def binningRequested = params.binning?.toString()?.toBoolean() ?: false
    def doBinning        = doAssembly && binningRequested
    if (!doAssembly && binningRequested) {
      log.warn "Warning: --binning is ignored because --noassembly was set"
    }

    // Only require assembly params when assembly is enabled
    if (doAssembly) {
      if (!params.taxdump || !params.uniprot_db) {
        log.error "Error: Missing --taxdump or --uniprot_db"
        missingParametersError()
      }
    }

    if (doScreening) {
      if (!params.taxdump || !params.gtdb_ncbi_map || !params.singlem_db) {
        log.error "Error: Missing --taxdump, --gtdb_ncbi_map, or --singlem_db required for taxa filtering"
        missingParametersError()
      }
      if (sraMode && !params.sandpiper_db) {
        log.error "Error: --sandpiper_db is required for SRA-based taxa screening"
        missingParametersError()
      }
    }
    else if (!doAssembly) {
      log.warn "Warning: --noassembly was set but --taxa was not provided; SINGLEM will not run in this mode."
    }


    // Common channels
    def outdir = file(params.outdir ?: './output').toAbsolutePath().toString()

    // taxdump is needed for screening and/or assembly; keep it optional when not needed
    def taxdump_ch = params.taxdump ? channel.value(file(params.taxdump)) : channel.empty()

    // uniprot is only needed for assembly (DIAMOND/SEMIBIN)
    def uniprot_db_ch = (doAssembly && params.uniprot_db) ? channel.value(file(params.uniprot_db)) : channel.empty()


    // Taxa-related channels
    def validated_taxa_ch = channel.empty()
    def singlem_db_ch     = channel.empty()
    def sandpiper_db_ch   = channel.empty()
    if (doScreening) {
      def taxa_ch          = channel.value( file(params.taxa) )
      def gtdb_ncbi_map_ch = channel.value( file(params.gtdb_ncbi_map) )
      singlem_db_ch        = channel.value( file(params.singlem_db) )
      sandpiper_db_ch      = sraMode ? channel.value( file(params.sandpiper_db) ) : channel.empty()

      // validate taxa
      validated_taxa_ch = VALIDATE_TAXA(taxa_ch, taxdump_ch, gtdb_ncbi_map_ch).valid_taxa
    }


    // SRA mode
    def sra_singlem_reads_ch    = channel.empty()
    def sra_metadata_skipped_ch = channel.empty()
    def sra_metadata_note       = channel.empty()
    def sra_sandpiper_note      = channel.empty()
    def sra_download_srr_note   = channel.empty()
    def sra_singlem_note        = channel.empty()
    if (sraMode) {
      sra_accessions_channel = channel.fromPath(params.sra, checkIfExists: true)
                      .splitCsv(header: true, strip: true)
                      .map { row -> row.sra.trim() }
                      .filter { row -> row }
                      .distinct()

      pre_screening_out = PRE_SCREENING(sra_accessions_channel, validated_taxa_ch, sandpiper_db_ch, singlem_db_ch)

      sra_singlem_reads_ch    = pre_screening_out.singlem_reads
      sra_metadata_skipped_ch = pre_screening_out.sra_metadata_skipped
      sra_metadata_note       = pre_screening_out.sra_metadata_note
      sra_sandpiper_note      = pre_screening_out.sandpiper_note
      sra_download_srr_note   = pre_screening_out.download_srr_note
      sra_singlem_note        = pre_screening_out.singlem_note
    }

    // FASTQ mode
    def fastq_singlem_reads_ch = channel.empty()
    def fastq_singlem_note     = channel.empty()
    if (fastqMode) {
      fastq_samplesheet_channel = channel.fromPath(params.fastq_tsv, checkIfExists: true)
                        .splitCsv(header: true, sep: '\t', strip: true)
                        .map { row ->
                          def sample    = (row.sample ?: '').trim()
                          def read_type = (row.read_type ?: '').trim()
                          def reads_raw = (row.reads ?: '').trim()

                          if (!sample || !read_type || !reads_raw) {
                            log.warn "Skipping FASTQ TSV row with missing fields: ${row}"
                            return null
                          }

                          def read_files = reads_raw.split(/\s*,\s*/).findAll { token -> token }.collect { read_path -> file(read_path) }

                          if (!read_files) {
                            log.warn "No valid FASTQ paths for sample ${sample}; skipping"
                            return null
                          }

                          def sra       = sample
                          def srr       = sample
                          def platform  = "UNKNOWN"
                          def model     = read_type
                          def strategy  = "UNKNOWN"

                          tuple(sra, srr, platform, model, strategy, read_type, read_files)
                        }
                        .filter { row -> row != null }

      def fastq_pre_screening_out = FASTQ_PRE_SCREENING(fastq_samplesheet_channel, validated_taxa_ch, singlem_db_ch)
      fastq_singlem_reads_ch = fastq_pre_screening_out.reads
      fastq_singlem_note     = fastq_pre_screening_out.note
    }

    // Merge SRA and FASTQ reads for assembly
    def screened_reads_channel = channel.empty()
                                  .mix(sra_singlem_reads_ch)
                                  .mix(fastq_singlem_reads_ch)
    def prescreening_note_channel = channel.empty()
                                  .mix(sra_singlem_note)
                                  .mix(fastq_singlem_note)


  // Empty notes
  def assembly_notes_ch = channel.empty()
  def diamond_note_ch   = channel.empty()
  def blobtools_note_ch = channel.empty()
  def taxa_note_ch      = channel.empty()
  def taxa_summary_ch   = channel.empty()
  def binning_note_entries_ch = channel.empty()

  if (noAssembly) {
    log.info "--noassembly set: skipping ASSEMBLY/BINNING; generating screening-only summary.tsv"

    // Prepare empty summary for succeeded samples
    def prescreen_success_meta = screened_reads_channel
      .map { sra, srr, platform, model, strategy, read_type, reads ->
        tuple(sra, srr, platform, model, strategy, read_type, '', '')
      }
      .distinct()
    def prescreen_empty = CREATE_EMPTY_SUMMARY(prescreen_success_meta).skipped_rows

    // Convert to the shape expected by SUMMARY's taxa_summary input:
    // (sra, srr, platform, model, strategy, read_type, assembler, summary_csv)
    taxa_summary_ch = prescreen_empty.map { sra, srr, platform, model, strategy, read_type, assembler, summary_csv, note ->
      tuple(sra, srr, platform, model, strategy, read_type, assembler, summary_csv)
    }
  }
  else {
    def assembly_out = ASSEMBLY(screened_reads_channel, validated_taxa_ch, uniprot_db_ch, taxdump_ch)

    assembly_notes_ch = assembly_out.assembly_notes
    diamond_note_ch   = assembly_out.diamond_note
    blobtools_note_ch = assembly_out.blobtools_note
    taxa_note_ch      = assembly_out.taxa_note
    taxa_summary_ch   = assembly_out.taxa_summary

    if (doBinning) {
      def binning_out = BINNING(assembly_out.blobtable, assembly_out.assembly_bam_all, uniprot_db_ch)
      binning_note_entries_ch = binning_out.note_entries
    }
  }


  SUMMARY(
    // PRE_SCREENING: summary-related outputs
    sra_metadata_skipped_ch,
    sra_metadata_note,
    sra_sandpiper_note,
    sra_download_srr_note,
    prescreening_note_channel,

    // ASSEMBLY: summary-related outputs (empty in --noassembly)
    assembly_notes_ch,
    diamond_note_ch,
    blobtools_note_ch,
    taxa_note_ch,
    taxa_summary_ch,

    // BINNING: raw note entries (empty in --noassembly or if binning off)
    binning_note_entries_ch,

    // constant
    outdir
  )

  workflow.onComplete = {
    def outdirPath  = file(params.outdir ?: './output').toAbsolutePath()
    def summaryFile = outdirPath.resolve('summary.tsv')
    def traceFile   = file("${workflow.launchDir}/execution-reports/trace.tsv").toAbsolutePath()
    def scriptFile  = file("${workflow.projectDir}/bin/annotate_summary_from_trace.py").toAbsolutePath()

    log.info "onComplete: summary.tsv -> ${summaryFile}"
    log.info "onComplete: trace.tsv   -> ${traceFile}"
    log.info "onComplete: annotator   -> ${scriptFile}"

    if( !summaryFile.exists() ) {
      log.warn "onComplete: ${summaryFile} not found; skipping scheduler annotation"
    }
    else if( !traceFile.exists() ) {
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
        .directory( workflow.launchDir.toFile() )
        .redirectError( java.lang.ProcessBuilder.Redirect.INHERIT )
        .redirectOutput( java.lang.ProcessBuilder.Redirect.INHERIT )
        .start()

      def rc = proc.waitFor()
      if( rc != 0 ) {
        log.warn "onComplete: annotator script exited with code ${rc}"
      }
      else {
        log.info "onComplete: summary.tsv successfully annotated with scheduler error information"
      }
    }
  }
}
