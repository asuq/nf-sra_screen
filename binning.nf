#!/usr/bin/env nextflow
nextflow.enable.dsl=2


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
    --max_retries   Maximum number of retries for each process (default: 3)
    --help          Show this help message

  binning.tsv columns:
    sample          Logical sample identifier
    read_type       short | nanopore | pacbio | hifi
    reads           Comma-separated FASTQ paths
    assembly_fasta  Path to assembly FASTA
  """.stripIndent()
}


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


/*
 * Parse and validate one binning samplesheet row.
 */
def normaliseBinningRow(row) {
  def sampleRaw = (row.sample ?: '').trim()
  def readTypeRaw = (row.read_type ?: '').trim().toLowerCase()
  def readsRaw = (row.reads ?: '').trim()
  def assemblyRaw = (row.assembly_fasta ?: '').trim()

  if (!sampleRaw || !readTypeRaw || !readsRaw || !assemblyRaw) {
    error "Invalid binning.tsv row with missing fields: ${row}"
  }

  def supportedReadTypes = ['short', 'nanopore', 'pacbio', 'hifi'] as Set
  if (!(readTypeRaw in supportedReadTypes)) {
    error "Unsupported read_type '${readTypeRaw}' for sample '${sampleRaw}'"
  }

  def readFiles = readsRaw
    .split(/\s*,\s*/)
    .findAll { it }
    .collect { file(it, checkIfExists: true) }

  if (!readFiles) {
    error "No read files were provided for sample '${sampleRaw}'"
  }

  if (readTypeRaw == 'short' && readFiles.size() > 2) {
    error "Sample '${sampleRaw}' has ${readFiles.size()} short-read FASTQs; expected one or two"
  }

  def assemblyFasta = file(assemblyRaw, checkIfExists: true)

  def sra = sampleRaw
  def srr = sampleRaw
  def platform = 'UNKNOWN'
  def model = readTypeRaw
  def strategy = 'UNKNOWN'
  def assembler = 'external_assembly'

  tuple(sra, srr, platform, model, strategy, assembler, assemblyFasta, readFiles)
}


process MAP_SHORT {
    tag "${sra}:${srr}"
    label 'binning'

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path(assembly_fasta), path(reads)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path(assembly_fasta), path("assembly.bam"), path("assembly.bam.csi"), optional: true, emit: mapped
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("FAIL.note"), optional: true, emit: note

    script:
    def mapScript = file("${workflow.projectDir}/bin/run_map_to_assembly.sh").toAbsolutePath()
    """
    ${mapScript} \\
      --read-type short \\
      --assembly "${assembly_fasta}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries} \\
      --reads ${reads}
    """

    stub:
    """
    : > assembly.bam
    : > assembly.bam.csi
    """
}


process MAP_NANO {
    tag "${sra}:${srr}"
    label 'binning'

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path(assembly_fasta), path(reads)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path(assembly_fasta), path("assembly.bam"), path("assembly.bam.csi"), optional: true, emit: mapped
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("FAIL.note"), optional: true, emit: note

    script:
    def mapScript = file("${workflow.projectDir}/bin/run_map_to_assembly.sh").toAbsolutePath()
    """
    ${mapScript} \\
      --read-type nanopore \\
      --assembly "${assembly_fasta}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries} \\
      --reads ${reads}
    """

    stub:
    """
    : > assembly.bam
    : > assembly.bam.csi
    """
}


process MAP_PACBIO {
    tag "${sra}:${srr}"
    label 'binning'

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path(assembly_fasta), path(reads)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path(assembly_fasta), path("assembly.bam"), path("assembly.bam.csi"), optional: true, emit: mapped
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("FAIL.note"), optional: true, emit: note

    script:
    def mapScript = file("${workflow.projectDir}/bin/run_map_to_assembly.sh").toAbsolutePath()
    """
    ${mapScript} \\
      --read-type pacbio \\
      --assembly "${assembly_fasta}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries} \\
      --reads ${reads}
    """

    stub:
    """
    : > assembly.bam
    : > assembly.bam.csi
    """
}


process MAP_HIFI {
    tag "${sra}:${srr}"
    label 'binning'

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path(assembly_fasta), path(reads)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path(assembly_fasta), path("assembly.bam"), path("assembly.bam.csi"), optional: true, emit: mapped
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("FAIL.note"), optional: true, emit: note

    script:
    def mapScript = file("${workflow.projectDir}/bin/run_map_to_assembly.sh").toAbsolutePath()
    """
    ${mapScript} \\
      --read-type hifi \\
      --assembly "${assembly_fasta}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries} \\
      --reads ${reads}
    """

    stub:
    """
    : > assembly.bam
    : > assembly.bam.csi
    """
}


process METABAT {
    tag "${sra}:${srr}"
    label 'binning'
    publishDir "${params.outdir}/${sra}/${srr}/binning",
      mode: 'copy',
      overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path(assembly_fasta), path(assembly_bam), path(assembly_csi)

    output:
    tuple val(sra), val(srr), path("metabat"),                                                                  emit: bins
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("metabat.note"), emit: note

    script:
    def metabatScript = file("${workflow.projectDir}/bin/run_metabat.sh").toAbsolutePath()
    """
    ${metabatScript} \\
      --assembly "${assembly_fasta}" \\
      --bam "${assembly_bam}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries}
    """

    stub:
    """
    mkdir -p metabat
    : > metabat.note
    """
}


process SEMIBIN {
    tag "${sra}:${srr}"
    label 'binning'
    publishDir "${params.outdir}/${sra}/${srr}/binning",
      mode: 'copy',
      overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path(assembly_fasta), path(assembly_bam), path(assembly_csi)
    path uniprot_db

    output:
    tuple val(sra), val(srr), path("semibin"), path("semibin/contig_bins.tsv"),                                 emit: bins
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("semibin.note"), emit: note

    script:
    def semibinScript = file("${workflow.projectDir}/bin/run_semibin.sh").toAbsolutePath()
    """
    ${semibinScript} \\
      --assembly "${assembly_fasta}" \\
      --bam "${assembly_bam}" \\
      --diamond-db "${uniprot_db}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries}
    """

    stub:
    """
    mkdir -p semibin
    printf 'contig\tbin\n' > semibin/contig_bins.tsv
    : > semibin.note
    """
}


process ROSELLA {
    tag "${sra}:${srr}"
    label 'binning'
    publishDir "${params.outdir}/${sra}/${srr}/binning",
      mode: 'copy',
      overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path(assembly_fasta), path(assembly_bam), path(assembly_csi)

    output:
    tuple val(sra), val(srr), path("rosella"),                                                                  emit: bins
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("rosella.note"), emit: note

    script:
    def rosellaScript = file("${workflow.projectDir}/bin/run_rosella.sh").toAbsolutePath()
    """
    ${rosellaScript} \\
      --assembly "${assembly_fasta}" \\
      --bam "${assembly_bam}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries}
    """

    stub:
    """
    mkdir -p rosella
    : > rosella.note
    """
}


process DASTOOL {
    tag "${sra}:${srr}"
    publishDir "${params.outdir}/${sra}/${srr}/binning",
      mode: 'copy',
      overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler),
          path(assembly_fasta),
          path(metabat_dir),
          path(semibin_dir),
          path(semibin_contig_bins),
          path(rosella_dir)

    output:
    tuple val(sra), val(srr), path("dastool"),                                                                  emit: bins
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("dastool.note"), emit: note

    script:
    def metaDir = metabat_dir ?: ''
    def semibinDir = semibin_dir ?: ''
    def semibinMap = semibin_contig_bins ?: ''
    def rosellaDir = rosella_dir ?: ''
    def dastoolScript = file("${workflow.projectDir}/bin/run_dastool.sh").toAbsolutePath()
    """
    ${dastoolScript} \\
      --assembly "${assembly_fasta}" \\
      --metabat-dir "${metaDir}" \\
      --semibin-dir "${semibinDir}" \\
      --semibin-map "${semibinMap}" \\
      --rosella-dir "${rosellaDir}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries}
    """

    stub:
    """
    mkdir -p dastool
    : > dastool.note
    """
}


process CREATE_EMPTY_SUMMARY {
    tag "${sra}:${srr}"

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), val(note)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("empty_summary.csv"), val(note), emit: skipped_rows

    script:
    """
    echo 'rank,ncbi_taxa,n_contigs,output_ids_csv,output_fasta' > empty_summary.csv
    """
}


process CREATE_EMPTY_FAILURE_SUMMARY {
    tag "${sra}:${srr}"

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), val(note)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("empty_summary.csv"), val(note), emit: skipped_rows

    script:
    """
    echo 'rank,ncbi_taxa,n_contigs,output_ids_csv,output_fasta' > empty_summary.csv
    """
}


process APPEND_SUMMARY {
    tag "${sra}:${srr}"

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path(summary_csv), val(note)
    val outdir

    output:
    path("summary.tsv"), optional: true, emit: global_summary

    script:
    def appendScript = file("${workflow.projectDir}/bin/run_append_summary.sh").toAbsolutePath()
    """
    ${appendScript} \\
      --outdir "${outdir}" \\
      --sra "${sra}" \\
      --srr "${srr}" \\
      --platform "${platform}" \\
      --model "${model}" \\
      --strategy "${strategy}" \\
      --assembler "${assembler}" \\
      --summary-csv "${summary_csv}" \\
      --note "${note}"
    """
}


workflow {
    if (params.help) {
      helpMessage()
      exit 0
    }

    if (!params.binning_tsv || !params.uniprot_db) {
      missingParametersError()
    }

    def outdir = file(params.outdir ?: './output').toAbsolutePath().toString()
    def uniprot_db_ch = channel.value(file(params.uniprot_db, checkIfExists: true))

    def binning_rows = channel.fromPath(params.binning_tsv, checkIfExists: true)
      .splitCsv(header: true, sep: '\t', strip: true)
      .map { row -> normaliseBinningRow(row) }

    def short_ch = binning_rows
      .filter { sra, srr, platform, model, strategy, assembler, assembly_fasta, reads ->
        model.equalsIgnoreCase('short')
      }
    def nano_ch = binning_rows
      .filter { sra, srr, platform, model, strategy, assembler, assembly_fasta, reads ->
        model.equalsIgnoreCase('nanopore')
      }
    def pacbio_ch = binning_rows
      .filter { sra, srr, platform, model, strategy, assembler, assembly_fasta, reads ->
        model.equalsIgnoreCase('pacbio')
      }
    def hifi_ch = binning_rows
      .filter { sra, srr, platform, model, strategy, assembler, assembly_fasta, reads ->
        model.equalsIgnoreCase('hifi')
      }

    def short_mapping = MAP_SHORT(short_ch)
    def nano_mapping = MAP_NANO(nano_ch)
    def pacbio_mapping = MAP_PACBIO(pacbio_ch)
    def hifi_mapping = MAP_HIFI(hifi_ch)

    def mapped_all = channel.empty()
      .mix(short_mapping.mapped)
      .mix(nano_mapping.mapped)
      .mix(pacbio_mapping.mapped)
      .mix(hifi_mapping.mapped)

    def mapping_note_ch = channel.empty()
      .mix(short_mapping.note)
      .mix(nano_mapping.note)
      .mix(pacbio_mapping.note)
      .mix(hifi_mapping.note)

    def metabat_binning = METABAT(mapped_all)
    def semibin_binning = SEMIBIN(mapped_all, uniprot_db_ch)
    def rosella_binning = ROSELLA(mapped_all)

    def dastool_base_by = mapped_all.map { sra, srr, platform, model, strategy, assembler, assembly_fasta, bam, csi ->
      tuple([sra, srr], [platform, model, strategy, assembler, assembly_fasta])
    }

    def metabat_entries = metabat_binning.bins.map { sra, srr, metabat_dir ->
      tuple([sra, srr], metabat_dir)
    }
    def semibin_entries = semibin_binning.bins.map { sra, srr, semibin_dir, semibin_bins ->
      tuple([sra, srr], [semibin_dir, semibin_bins])
    }
    def rosella_entries = rosella_binning.bins.map { sra, srr, rosella_dir ->
      tuple([sra, srr], rosella_dir)
    }

    def dastool_join = dastool_base_by
      .join(metabat_entries)
      .join(semibin_entries)
      .join(rosella_entries)

    def dastool_in = dastool_join.map { key, meta, metabat_dir, semibin_pair, rosella_dir ->
      def (sra, srr) = key
      def (platform, model, strategy, assembler, assembly_fasta) = meta
      def (semibin_dir, semibin_bins) = semibin_pair
      tuple(
        sra, srr, platform, model, strategy, assembler, assembly_fasta,
        metabat_dir, semibin_dir, semibin_bins, rosella_dir
      )
    }

    def dastool_binning = DASTOOL(dastool_in)

    def success_meta = mapped_all
      .map { sra, srr, platform, model, strategy, assembler, assembly_fasta, bam, csi ->
        tuple(sra, srr, platform, model, strategy, assembler, '')
      }
      .distinct()

    def empty_summary = CREATE_EMPTY_SUMMARY(success_meta).skipped_rows
    def succeeded_rows = empty_summary.map { sra, srr, platform, model, strategy, assembler, summary_csv, note ->
      tuple(sra, srr, platform, model, strategy, assembler, summary_csv, '')
    }

    def binning_notes = channel.empty()
      .mix(metabat_binning.note.map { sra, srr, platform, model, strategy, assembler, note_path ->
        tuple([sra, srr, platform, model, strategy, assembler], ['metabat', note_path])
      })
      .mix(semibin_binning.note.map { sra, srr, platform, model, strategy, assembler, note_path ->
        tuple([sra, srr, platform, model, strategy, assembler], ['semibin', note_path])
      })
      .mix(rosella_binning.note.map { sra, srr, platform, model, strategy, assembler, note_path ->
        tuple([sra, srr, platform, model, strategy, assembler], ['rosella', note_path])
      })
      .mix(dastool_binning.note.map { sra, srr, platform, model, strategy, assembler, note_path ->
        tuple([sra, srr, platform, model, strategy, assembler], ['dastool', note_path])
      })
      .groupTuple()
      .map { key, items ->
        def (sra, srr, platform, model, strategy, assembler) = key
        def noteText = items
          .sort { left, right -> left[0] <=> right[0] }
          .collect { entry ->
            def notePath = entry[1]
            def text = file(notePath as String).text.trim()
            text ? "${entry[0]}: ${text}" : null
          }
          .findAll { it }
          .join('; ')
        tuple(sra, srr, platform, model, strategy, assembler, noteText)
      }

    def final_success = succeeded_rows
      .map { sra, srr, platform, model, strategy, assembler, summary_csv, note ->
        tuple([sra, srr, platform, model, strategy, assembler], [summary_csv, note])
      }
      .join(binning_notes.map { sra, srr, platform, model, strategy, assembler, note ->
        tuple([sra, srr, platform, model, strategy, assembler], note)
      })
      .map { key, summary_entry, binning_note ->
        def (sra, srr, platform, model, strategy, assembler) = key
        def summary_csv = summary_entry[0]
        def base_note = summary_entry[1] ?: ''
        def extra_note = (binning_note ?: '').trim()
        def final_note = base_note
        if (extra_note) {
          final_note = final_note ? "${final_note}; ${extra_note}" : extra_note
        }
        tuple(sra, srr, platform, model, strategy, assembler, summary_csv, final_note)
      }

    def mapping_errors = mapping_note_ch.map { sra, srr, platform, model, strategy, assembler, note_path ->
      def note = file(note_path).text.trim()
      tuple(sra, srr, platform, model, strategy, assembler, note)
    }

    def failed_rows = CREATE_EMPTY_FAILURE_SUMMARY(mapping_errors)

    def summary_rows = channel.empty()
      .mix(final_success)
      .mix(failed_rows)

    APPEND_SUMMARY(summary_rows, outdir)
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
        return
      }

      if (!traceFile.exists()) {
        log.warn "onComplete: ${traceFile} not found; skipping scheduler annotation"
        return
      }

      def cmd = [
        'python3',
        scriptFile.toString(),
        summaryFile.toString(),
        traceFile.toString()
      ]

      log.info "onComplete: running ${cmd.join(' ')}"

      def proc = new ProcessBuilder(cmd)
        .directory(workflow.launchDir.toFile())
        .redirectError(java.lang.ProcessBuilder.Redirect.INHERIT)
        .redirectOutput(java.lang.ProcessBuilder.Redirect.INHERIT)
        .start()

      int rc = proc.waitFor()
      if (rc != 0) {
        log.warn "onComplete: annotator script exited with code ${rc}"
      }
      else {
        log.info "onComplete: summary.tsv successfully annotated with scheduler error information"
      }
    }
}
