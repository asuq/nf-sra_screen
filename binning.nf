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
    --binners       Comma-separated binners (default: auto; allowed: auto,metabat,semibin,rosella,comebin,vamb,lorbin)
    --refiners      Comma-separated refiners (default: dastool; allowed: dastool,binette)
    --checkm2_db    CheckM2 DIAMOND database required with --refiners binette
    --semibin_environment  SemiBin2 pretrained environment (default: global)
    --gpu           Use GPU variants for COMEBin, VAMB, and HiFi-only LorBin
    --gpu_type      GPU type for scheduler requests on GWDG (default: A100)
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


/*
 * Return all implemented binners in their stable execution order.
 */
def allImplementedBinners() {
  return ['metabat', 'semibin', 'rosella', 'comebin', 'vamb', 'lorbin']
}


/*
 * Return the binners compatible with the supplied read type.
 */
def compatibleBinnersForReadType(readType) {
  def binners = ['metabat', 'semibin', 'rosella', 'comebin', 'vamb']
  if (readType?.toString()?.trim()?.equalsIgnoreCase('hifi')) {
    binners = binners + ['lorbin']
  }
  return binners
}


/*
 * Parse and validate a comma-separated binning tool selection.
 */
def parseBinnerSelection(rawValue) {
  def raw = rawValue == null ? 'auto' : rawValue.toString().trim().toLowerCase()
  if (!raw) {
    error "--binners must include at least one tool"
  }

  if (raw in ['auto', 'all-compatible']) {
    return [mode: 'auto', tools: allImplementedBinners()]
  }

  return [
    mode: 'explicit',
    tools: parsePhase0ToolSelection(
      raw,
      'auto',
      allImplementedBinners() as Set,
      [] as Set,
      'binners'
    )
  ]
}


/*
 * Return the effective binners for one sample, or fail on invalid LorBin use.
 */
def binnersForSample(selection, sra, srr, readType) {
  def readTypeLc = readType?.toString()?.trim()?.toLowerCase()
  def mode = selection.mode ?: selection.binnerMode
  def tools = selection.tools ?: selection.binners
  def selected = mode == 'auto'
    ? compatibleBinnersForReadType(readTypeLc)
    : tools

  if ('lorbin' in selected && readTypeLc != 'hifi') {
    error "LorBin only supports hifi reads; sample ${sra}:${srr} has read_type ${readType}. Remove lorbin from --binners for this run."
  }

  return selected
}


/*
 * Test whether a comma-separated binner list includes one tool.
 */
def binnerCsvContains(csv, tool) {
  return csv.toString().split(',').collect { it.trim() }.contains(tool)
}


/*
 * Count tools in a comma-separated binner list.
 */
def binnerCsvSize(csv) {
  return csv.toString().split(',').collect { it.trim() }.findAll { it }.size()
}


/*
 * Parse and validate a comma-separated tool selection.
 */
def parsePhase0ToolSelection(rawValue, defaultValue, allowedTools, plannedTools, paramName) {
  def raw = rawValue == null ? defaultValue : rawValue.toString()
  def selected = raw
    .split(',')
    .collect { it.trim().toLowerCase() }
    .findAll { it }
    .unique()

  if (!selected) {
    error "--${paramName} must include at least one tool"
  }

  def planned = selected.findAll { it in plannedTools }
  if (planned) {
    error "--${paramName} includes planned tool(s) not implemented yet: ${planned.join(', ')}"
  }

  def invalid = selected.findAll { !(it in allowedTools) }
  if (invalid) {
    error "--${paramName} includes unsupported tool(s): ${invalid.join(', ')}"
  }

  selected
}


/*
 * Validate the selected SemiBin2 pretrained environment.
 */
def validateSemibinEnvironment(rawValue) {
  def selected = (rawValue ?: 'global').toString().trim().toLowerCase()
  def allowed = [
    'human_gut',
    'dog_gut',
    'ocean',
    'soil',
    'cat_gut',
    'human_oral',
    'mouse_gut',
    'pig_gut',
    'built_environment',
    'wastewater',
    'chicken_caecum',
    'global'
  ] as Set

  if (!selected) {
    error "--semibin_environment must not be empty"
  }

  if (!(selected in allowed)) {
    error "--semibin_environment includes unsupported environment: ${selected}"
  }

  selected
}


/*
 * Validate binning syntax before workflow construction.
 */
def validatePhase0BinningOptions() {
  def useGpu = params.gpu?.toString()?.toBoolean() ?: false
  def plannedTools = [] as Set
  def binnerSelection = parseBinnerSelection(params.binners)
  def refiners = parsePhase0ToolSelection(
    params.refiners,
    'dastool',
    ['dastool', 'binette'] as Set,
    plannedTools,
    'refiners'
  )

  if ('binette' in refiners && !params.checkm2_db) {
    error "--checkm2_db is required when --refiners includes binette"
  }

  [
    binners: binnerSelection.tools,
    binnerMode: binnerSelection.mode,
    refiners: refiners,
    gpu: useGpu,
    semibinEnvironment: validateSemibinEnvironment(params.semibin_environment)
  ]
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
  def sampleRaw = (row['sample'] ?: '').trim()
  def srrRaw = (row['srr'] ?: '').trim()
  def readTypeRaw = (row['read_type'] ?: '').trim().toLowerCase()
  def readsRaw = (row['reads'] ?: '').trim()
  def assemblyRaw = (row['assembly_fasta'] ?: '').trim()

  if (!sampleRaw || !assemblyRaw) {
    error "Invalid binning.tsv row with missing sample or assembly_fasta: ${row}"
  }

  def hasReads = !!readsRaw
  def hasSrr = !!srrRaw
  if (hasReads == hasSrr) {
    error "Sample '${sampleRaw}' must provide exactly one of reads or srr"
  }

  def assemblyFasta = file(assemblyRaw, checkIfExists: true)

  if (hasSrr) {
    return [
      mode: 'srr',
      sra: sampleRaw,
      srr: srrRaw,
      assembly_fasta: assemblyFasta
    ]
  }

  if (!readTypeRaw) {
    error "Sample '${sampleRaw}' is missing read_type for local reads"
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

  return [
    mode: 'local',
    sra: sampleRaw,
    srr: sampleRaw,
    platform: 'UNKNOWN',
    model: 'UNKNOWN',
    strategy: 'UNKNOWN',
    read_type: readTypeRaw,
    assembler: 'provided',
    assembly_fasta: assemblyFasta,
    reads: readFiles
  ]
}


/*
 * Parse the resolved SRR metadata emitted by the resolver helper.
 */
def parseResolvedMetadata(sra, srr, resolvedTsv) {
  def fields = file(resolvedTsv).text
    .readLines()
    .find { it?.trim() }
    ?.split(/\t/, -1)
    ?.collect { it.trim() }

  if (!fields || fields.size() != 4) {
    error "Invalid resolved metadata for sample '${sra}' run '${srr}'"
  }

  def (platform, model, strategy, readType) = fields
  if (!platform || !model || !strategy || !readType) {
    error "Incomplete resolved metadata for sample '${sra}' run '${srr}'"
  }

  tuple(sra, srr, platform, model, strategy, readType)
}


process RESOLVE_SRR_METADATA {
    tag "${sra}:${srr}"

    input:
    tuple val(sra), val(srr), path(assembly_fasta)

    output:
    tuple val(sra), val(srr), path("resolved.tsv"), path(assembly_fasta), optional: true, emit: resolved
    tuple val(sra), val(srr), val('UNKNOWN'), val('UNKNOWN'), val('UNKNOWN'), val('UNKNOWN'), val('provided'), path("FAIL.note"), optional: true, emit: note

    script:
    def resolveScript = file("${workflow.projectDir}/bin/run_resolve_srr_metadata.sh").toAbsolutePath()
    """
    ${resolveScript} \\
      --sample "${sra}" \\
      --srr "${srr}" \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries}
    """

    stub:
    """
    printf 'ILLUMINA\tNovaSeq 6000\tWGS\tshort\n' > resolved.tsv
    """
}


process DOWNLOAD_SRR {
    tag "${sra}:${srr}"
    publishDir "${params.outdir}/${sra}/${srr}/",
      mode: 'copy',
      overwrite: true,
      saveAs: { filename ->
        filename in ["FAIL.note"] ? filename : null
      }

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(assembly_fasta)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(assembly_fasta), path("*.f*q*"), path("assembler.txt"), optional: true, emit: reads
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("FAIL.note"), optional: true, emit: note

    script:
    """
    run_download_srr.sh \\
      --srr "${srr}" \\
      --platform "${platform}" \\
      --read-type "${read_type}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries}
    """

    stub:
    """
    : > "${srr}_1.fastq.gz"
    printf '%s\n' "${read_type}" > assembler.txt
    """
}


process MAP_SHORT {
    tag "${sra}:${srr}"
    label 'binning'

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(assembly_fasta), path(reads)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("assembly.bam"), path("assembly.bam.csi"), optional: true, emit: mapped
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("FAIL.note"), optional: true, emit: note

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
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(assembly_fasta), path(reads)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("assembly.bam"), path("assembly.bam.csi"), optional: true, emit: mapped
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("FAIL.note"), optional: true, emit: note

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
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(assembly_fasta), path(reads)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("assembly.bam"), path("assembly.bam.csi"), optional: true, emit: mapped
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("FAIL.note"), optional: true, emit: note

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
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(assembly_fasta), path(reads)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("assembly.bam"), path("assembly.bam.csi"), optional: true, emit: mapped
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("FAIL.note"), optional: true, emit: note

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
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(assembly_fasta), path(assembly_bam), path(assembly_csi)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler),
          val("metabat"), path("metabat"), path("metabat.contig2bin.tsv"), path("metabat.note"),               emit: result

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
    : > metabat.contig2bin.tsv
    : > metabat.note
    """
}


process COMEBIN {
    tag "${sra}:${srr}"
    label 'binning'
    publishDir "${params.outdir}/${sra}/${srr}/binning",
      mode: 'copy',
      overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(assembly_fasta), path(assembly_bam), path(assembly_csi)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler),
          val("comebin"), path("comebin"), path("comebin.contig2bin.tsv"), path("comebin.note"),             emit: result

    script:
    def comebinScript = file("${workflow.projectDir}/bin/run_comebin_nf.sh").toAbsolutePath()
    """
    ${comebinScript} \\
      --assembly "${assembly_fasta}" \\
      --bam "${assembly_bam}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries}
    """

    stub:
    """
    mkdir -p comebin
    : > comebin.contig2bin.tsv
    : > comebin.note
    """
}


process COMEBIN_GPU {
    tag "${sra}:${srr}"
    label 'gpu'
    publishDir "${params.outdir}/${sra}/${srr}/binning",
      mode: 'copy',
      overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(assembly_fasta), path(assembly_bam), path(assembly_csi)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler),
          val("comebin"), path("comebin"), path("comebin.contig2bin.tsv"), path("comebin.note"),             emit: result

    script:
    def comebinScript = file("${workflow.projectDir}/bin/run_comebin_nf.sh").toAbsolutePath()
    """
    ${comebinScript} \\
      --assembly "${assembly_fasta}" \\
      --bam "${assembly_bam}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries} \\
      --require-cuda
    """

    stub:
    """
    mkdir -p comebin
    : > comebin.contig2bin.tsv
    : > comebin.note
    """
}


process VAMB {
    tag "${sra}:${srr}"
    label 'binning'
    publishDir "${params.outdir}/${sra}/${srr}/binning",
      mode: 'copy',
      overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(assembly_fasta), path(assembly_bam), path(assembly_csi)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler),
          val("vamb"), path("vamb"), path("vamb.contig2bin.tsv"), path("vamb.note"),                         emit: result

    script:
    def vambScript = file("${workflow.projectDir}/bin/run_vamb.sh").toAbsolutePath()
    """
    ${vambScript} \\
      --assembly "${assembly_fasta}" \\
      --bam "${assembly_bam}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries}
    """

    stub:
    """
    mkdir -p vamb
    : > vamb.contig2bin.tsv
    : > vamb.note
    """
}


process VAMB_GPU {
    tag "${sra}:${srr}"
    label 'gpu'
    publishDir "${params.outdir}/${sra}/${srr}/binning",
      mode: 'copy',
      overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(assembly_fasta), path(assembly_bam), path(assembly_csi)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler),
          val("vamb"), path("vamb"), path("vamb.contig2bin.tsv"), path("vamb.note"),                         emit: result

    script:
    def vambScript = file("${workflow.projectDir}/bin/run_vamb.sh").toAbsolutePath()
    """
    ${vambScript} \\
      --assembly "${assembly_fasta}" \\
      --bam "${assembly_bam}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries} \\
      --cuda \\
      --require-cuda
    """

    stub:
    """
    mkdir -p vamb
    : > vamb.contig2bin.tsv
    : > vamb.note
    """
}


process LORBIN {
    tag "${sra}:${srr}"
    label 'binning'
    publishDir "${params.outdir}/${sra}/${srr}/binning",
      mode: 'copy',
      overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(assembly_fasta), path(assembly_bam), path(assembly_csi)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler),
          val("lorbin"), path("lorbin"), path("lorbin.contig2bin.tsv"), path("lorbin.note"),                  emit: result

    script:
    def lorbinScript = file("${workflow.projectDir}/bin/run_lorbin.sh").toAbsolutePath()
    """
    ${lorbinScript} \\
      --assembly "${assembly_fasta}" \\
      --bam "${assembly_bam}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries}
    """

    stub:
    """
    mkdir -p lorbin
    : > lorbin.contig2bin.tsv
    : > lorbin.note
    """
}


process LORBIN_GPU {
    tag "${sra}:${srr}"
    label 'gpu'
    publishDir "${params.outdir}/${sra}/${srr}/binning",
      mode: 'copy',
      overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(assembly_fasta), path(assembly_bam), path(assembly_csi)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler),
          val("lorbin"), path("lorbin"), path("lorbin.contig2bin.tsv"), path("lorbin.note"),                  emit: result

    script:
    def lorbinScript = file("${workflow.projectDir}/bin/run_lorbin.sh").toAbsolutePath()
    """
    ${lorbinScript} \\
      --assembly "${assembly_fasta}" \\
      --bam "${assembly_bam}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries} \\
      --require-cuda
    """

    stub:
    """
    mkdir -p lorbin
    : > lorbin.contig2bin.tsv
    : > lorbin.note
    """
}


process SEMIBIN {
    tag "${sra}:${srr}"
    label 'binning'
    publishDir "${params.outdir}/${sra}/${srr}/binning",
      mode: 'copy',
      overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(assembly_fasta), path(assembly_bam), path(assembly_csi)
    path uniprot_db

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler),
          val("semibin"), path("semibin"), path("semibin.contig2bin.tsv"), path("semibin.note"),               emit: result

    script:
    def semibinScript = file("${workflow.projectDir}/bin/run_semibin.sh").toAbsolutePath()
    """
    ${semibinScript} \\
      --assembly "${assembly_fasta}" \\
      --bam "${assembly_bam}" \\
      --diamond-db "${uniprot_db}" \\
      --read-type "${read_type}" \\
      --environment "${params.semibin_environment}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries}
    """

    stub:
    """
    mkdir -p semibin
    printf 'contig\tbin\n' > semibin/contig_bins.tsv
    : > semibin.contig2bin.tsv
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
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(assembly_fasta), path(assembly_bam), path(assembly_csi)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler),
          val("rosella"), path("rosella"), path("rosella.contig2bin.tsv"), path("rosella.note"),               emit: result

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
    : > rosella.contig2bin.tsv
    : > rosella.note
    """
}


process DASTOOL {
    tag "${sra}:${srr}"
    publishDir "${params.outdir}/${sra}/${srr}/binning",
      mode: 'copy',
      overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler),
          path(assembly_fasta),
          path(contig2bin_maps)

    output:
    tuple val(sra), val(srr), val(read_type), val(assembler), path("dastool"),                                                 emit: bins
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("dastool.note"), emit: note

    script:
    def mapArg = contig2bin_maps.collect { it.toString() }.join(',')
    def dastoolScript = file("${workflow.projectDir}/bin/run_dastool.sh").toAbsolutePath()
    """
    ${dastoolScript} \\
      --assembly "${assembly_fasta}" \\
      --bin-maps "${mapArg}" \\
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


process BINETTE {
    tag "${sra}:${srr}"
    publishDir "${params.outdir}/${sra}/${srr}/binning",
      mode: 'copy',
      overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler),
          path(assembly_fasta),
          path(contig2bin_maps)
    path checkm2_db

    output:
    tuple val(sra), val(srr), val(read_type), val(assembler), path("binette"),                                                 emit: bins
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("binette.note"), emit: note

    script:
    def mapArg = contig2bin_maps.collect { it.toString() }.join(',')
    def binetteScript = file("${workflow.projectDir}/bin/run_binette.sh").toAbsolutePath()
    """
    ${binetteScript} \\
      --assembly "${assembly_fasta}" \\
      --bin-maps "${mapArg}" \\
      --checkm2-db "${checkm2_db}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries}
    """

    stub:
    """
    mkdir -p binette
    : > binette/final_contig_to_bin.tsv
    : > binette.note
    """
}


process CREATE_EMPTY_SUMMARY {
    tag "${sra}:${srr}"

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), val(note)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("empty_summary.csv"), val(note), emit: skipped_rows

    script:
    """
    echo 'rank,ncbi_taxa,n_contigs,output_ids_csv,output_fasta' > empty_summary.csv
    """
}


process CREATE_EMPTY_FAILURE_SUMMARY {
    tag "${sra}:${srr}"

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), val(note)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("empty_summary.csv"), val(note), emit: skipped_rows

    script:
    """
    echo 'rank,ncbi_taxa,n_contigs,output_ids_csv,output_fasta' > empty_summary.csv
    """
}


process APPEND_SUMMARY {
    tag "${sra}:${srr}"

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(summary_csv), val(note)
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
      --read-type "${read_type}" \\
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

    def phase0Options = validatePhase0BinningOptions()
    def selectedBinners = phase0Options.binners
    def selectedRefiners = phase0Options.refiners
    def useGpu = phase0Options.gpu

    if (useGpu) {
      def cpuOnlyBinners = selectedBinners.findAll { it in ['metabat', 'rosella', 'semibin'] }
      if (cpuOnlyBinners) {
        log.warn "GPU mode requested; CPU-only binner(s) will remain on CPU: ${cpuOnlyBinners.join(', ')}"
      }
    }

    def outdir = file(params.outdir ?: './output').toAbsolutePath().toString()
    def uniprot_db_ch = channel.value(file(params.uniprot_db, checkIfExists: true))

    def binning_rows = channel.fromPath(params.binning_tsv, checkIfExists: true)
      .splitCsv(header: true, sep: '\t', strip: true)
      .map { row -> normaliseBinningRow(row) }

    def local_rows = binning_rows
      .filter { row -> row.mode == 'local' }
      .map { row ->
        tuple(
          row.sra,
          row.srr,
          row.platform,
          row.model,
          row.strategy,
          row.read_type,
          row.assembler,
          row.assembly_fasta,
          row.reads
        )
      }

    def srr_rows = binning_rows
      .filter { row -> row.mode == 'srr' }
      .map { row ->
        tuple(row.sra, row.srr, row.assembly_fasta)
      }

    def resolved_srr = RESOLVE_SRR_METADATA(srr_rows)
    def resolved_srr_rows = resolved_srr.resolved.map { sra, srr, resolved_tsv, assembly_fasta ->
      def resolved = parseResolvedMetadata(sra, srr, resolved_tsv)
      tuple(
        resolved[0],
        resolved[1],
        resolved[2],
        resolved[3],
        resolved[4],
        resolved[5],
        'provided',
        assembly_fasta
      )
    }

    def downloaded_srr = DOWNLOAD_SRR(resolved_srr_rows)
    def downloaded_rows = downloaded_srr.reads.map { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, reads, asm_txt ->
      def fixedReadType = file(asm_txt).text.trim() ?: read_type
      tuple(sra, srr, platform, model, strategy, fixedReadType, assembler, assembly_fasta, reads)
    }

    def mapping_rows = channel.empty()
      .mix(local_rows)
      .mix(downloaded_rows)

    def short_ch = mapping_rows
      .filter { _sra, _srr, _platform, _model, _strategy, read_type, _assembler, _assembly_fasta, _reads ->
        read_type.equalsIgnoreCase('short')
      }
    def nano_ch = mapping_rows
      .filter { _sra, _srr, _platform, _model, _strategy, read_type, _assembler, _assembly_fasta, _reads ->
        read_type.equalsIgnoreCase('nanopore')
      }
    def pacbio_ch = mapping_rows
      .filter { _sra, _srr, _platform, _model, _strategy, read_type, _assembler, _assembly_fasta, _reads ->
        read_type.equalsIgnoreCase('pacbio')
      }
    def hifi_ch = mapping_rows
      .filter { _sra, _srr, _platform, _model, _strategy, read_type, _assembler, _assembly_fasta, _reads ->
        read_type.equalsIgnoreCase('hifi')
      }

    def short_mapping = MAP_SHORT(short_ch)
    def nano_mapping = MAP_NANO(nano_ch)
    def pacbio_mapping = MAP_PACBIO(pacbio_ch)
    def hifi_mapping = MAP_HIFI(hifi_ch)

    def mapped_bams = channel.empty()
      .mix(short_mapping.mapped)
      .mix(nano_mapping.mapped)
      .mix(pacbio_mapping.mapped)
      .mix(hifi_mapping.mapped)

    def assembly_by_sample = mapping_rows.map { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, _reads ->
      tuple([sra, srr, read_type, assembler], [platform, model, strategy, assembly_fasta])
    }

    def mapped_bams_by_sample = mapped_bams.map { sra, srr, _platform, _model, _strategy, read_type, assembler, assembly_bam, assembly_csi ->
      tuple([sra, srr, read_type, assembler], [assembly_bam, assembly_csi])
    }

    def mapped_all = assembly_by_sample
      .join(mapped_bams_by_sample)
      .map { key, meta, bam_idx ->
        def (sra, srr, read_type, assembler) = key
        def (platform, model, strategy, assembly_fasta) = meta
        def (assembly_bam, assembly_csi) = bam_idx
        tuple(sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi)
      }

    def failure_note_ch = channel.empty()
      .mix(resolved_srr.note)
      .mix(downloaded_srr.note)
      .mix(short_mapping.note)
      .mix(nano_mapping.note)
      .mix(pacbio_mapping.note)
      .mix(hifi_mapping.note)

    def planned_binning_input = mapped_all
      .map { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi ->
        def sampleBinners = binnersForSample(phase0Options, sra, srr, read_type).join(',')
        tuple(sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi, sampleBinners)
      }

    def binner_plan_by_sample = planned_binning_input
      .map { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi, sampleBinners ->
        tuple([sra, srr, read_type, assembler], binnerCsvSize(sampleBinners))
      }

    def metabat_input = planned_binning_input
      .filter { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi, sampleBinners ->
        binnerCsvContains(sampleBinners, 'metabat')
      }
      .map { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi, sampleBinners ->
        tuple(sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi)
      }

    def comebin_input = planned_binning_input
      .filter { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi, sampleBinners ->
        binnerCsvContains(sampleBinners, 'comebin')
      }
      .map { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi, sampleBinners ->
        tuple(sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi)
      }

    def vamb_input = planned_binning_input
      .filter { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi, sampleBinners ->
        binnerCsvContains(sampleBinners, 'vamb')
      }
      .map { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi, sampleBinners ->
        tuple(sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi)
      }

    def lorbin_input = planned_binning_input
      .filter { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi, sampleBinners ->
        binnerCsvContains(sampleBinners, 'lorbin')
      }
      .map { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi, sampleBinners ->
        tuple(sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi)
      }

    def semibin_input = planned_binning_input
      .filter { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi, sampleBinners ->
        binnerCsvContains(sampleBinners, 'semibin')
      }
      .map { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi, sampleBinners ->
        tuple(sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi)
      }

    def rosella_input = planned_binning_input
      .filter { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi, sampleBinners ->
        binnerCsvContains(sampleBinners, 'rosella')
      }
      .map { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi, sampleBinners ->
        tuple(sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi)
      }

    def metabat_results = channel.empty()
    def comebin_results = channel.empty()
    def vamb_results = channel.empty()
    def lorbin_results = channel.empty()
    def semibin_results = channel.empty()
    def rosella_results = channel.empty()

    if ('metabat' in selectedBinners) {
      metabat_results = METABAT(metabat_input).result
    }
    if ('comebin' in selectedBinners) {
      if (useGpu) {
        comebin_results = COMEBIN_GPU(comebin_input).result
      }
      else {
        comebin_results = COMEBIN(comebin_input).result
      }
    }
    if ('vamb' in selectedBinners) {
      if (useGpu) {
        vamb_results = VAMB_GPU(vamb_input).result
      }
      else {
        vamb_results = VAMB(vamb_input).result
      }
    }
    if ('lorbin' in selectedBinners) {
      if (useGpu) {
        lorbin_results = LORBIN_GPU(lorbin_input).result
      }
      else {
        lorbin_results = LORBIN(lorbin_input).result
      }
    }
    if ('semibin' in selectedBinners) {
      semibin_results = SEMIBIN(semibin_input, uniprot_db_ch).result
    }
    if ('rosella' in selectedBinners) {
      rosella_results = ROSELLA(rosella_input).result
    }

    def binner_results = channel.empty()
      .mix(metabat_results)
      .mix(comebin_results)
      .mix(vamb_results)
      .mix(lorbin_results)
      .mix(semibin_results)
      .mix(rosella_results)

    def dastool_base_by = mapped_all.map { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, _bam, _csi ->
      tuple([sra, srr, read_type, assembler], [platform, model, strategy, assembly_fasta])
    }

    def binner_maps_by_sample = binner_results
      .map { sra, srr, platform, model, strategy, read_type, assembler, tool, bin_dir, contig2bin, note_path ->
        tuple([sra, srr, read_type, assembler], [tool, contig2bin])
      }
      .combine(binner_plan_by_sample, by: 0)
      .map { key, entry, expectedCount -> tuple(groupKey(key, expectedCount as int), entry) }
      .groupTuple()
      .map { key, entries -> tuple(key.target, entries) }

    def dastool_join = dastool_base_by.join(binner_maps_by_sample)

    def dastool_in = dastool_join.map { key, meta, entries ->
      def (sra, srr, read_type, assembler) = key
      def (platform, model, strategy, assembly_fasta) = meta
      def contig2bin_maps = entries.collect { entry -> entry[1] }
      tuple(
        sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta,
        contig2bin_maps
      )
    }

    def dastool_note_entries = channel.empty()
    if ('dastool' in selectedRefiners) {
      dastool_note_entries = DASTOOL(dastool_in).note.map { sra, srr, platform, model, strategy, read_type, assembler, note_path ->
        tuple(sra, srr, platform, model, strategy, read_type, assembler, 'dastool', note_path)
      }
    }

    def binette_note_entries = channel.empty()
    if ('binette' in selectedRefiners) {
      def checkm2_db_ch = Channel.value(file(params.checkm2_db, checkIfExists: true))
      binette_note_entries = BINETTE(dastool_in, checkm2_db_ch).note.map { sra, srr, platform, model, strategy, read_type, assembler, note_path ->
        tuple(sra, srr, platform, model, strategy, read_type, assembler, 'binette', note_path)
      }
    }

    def success_meta = mapped_all
      .map { sra, srr, platform, model, strategy, read_type, assembler, _assembly_fasta, _bam, _csi ->
        tuple(sra, srr, platform, model, strategy, read_type, assembler, '')
      }
      .distinct()

    def empty_summary = CREATE_EMPTY_SUMMARY(success_meta).skipped_rows
    def succeeded_rows = empty_summary.map { sra, srr, platform, model, strategy, read_type, assembler, summary_csv, _note ->
      tuple(sra, srr, platform, model, strategy, read_type, assembler, summary_csv, '')
    }

    def binner_note_entries = binner_results.map { sra, srr, platform, model, strategy, read_type, assembler, tool, bin_dir, contig2bin, note_path ->
      tuple(sra, srr, platform, model, strategy, read_type, assembler, tool, note_path)
    }

    def binning_notes = channel.empty()
      .mix(binner_note_entries)
      .mix(dastool_note_entries)
      .mix(binette_note_entries)
      .map { sra, srr, platform, model, strategy, read_type, assembler, tool, note_path ->
        tuple([sra, srr, platform, model, strategy, read_type, assembler], [tool, note_path])
      }
      .groupTuple()
      .map { key, items ->
        def (sra, srr, platform, model, strategy, read_type, assembler) = key
        def noteText = items
          .sort { left, right -> left[0] <=> right[0] }
          .collect { entry ->
            def notePath = entry[1]
            def text = file(notePath as String).text.trim()
            text ? "${entry[0]}: ${text}" : null
          }
          .findAll { it }
          .join('; ')
        tuple(sra, srr, platform, model, strategy, read_type, assembler, noteText)
      }

    def final_success = succeeded_rows
      .map { sra, srr, platform, model, strategy, read_type, assembler, summary_csv, note ->
        tuple([sra, srr, platform, model, strategy, read_type, assembler], [summary_csv, note])
      }
      .join(binning_notes.map { sra, srr, platform, model, strategy, read_type, assembler, note ->
        tuple([sra, srr, platform, model, strategy, read_type, assembler], note)
      })
      .map { key, summary_entry, binning_note ->
        def (sra, srr, platform, model, strategy, read_type, assembler) = key
        def summary_csv = summary_entry[0]
        def base_note = summary_entry[1] ?: ''
        def extra_note = (binning_note ?: '').trim()
        def final_note = base_note
        if (extra_note) {
          final_note = final_note ? "${final_note}; ${extra_note}" : extra_note
        }
        tuple(sra, srr, platform, model, strategy, read_type, assembler, summary_csv, final_note)
      }

    def mapping_errors = failure_note_ch.map { sra, srr, platform, model, strategy, read_type, assembler, note_path ->
      def note = file(note_path).text.trim()
      tuple(sra, srr, platform, model, strategy, read_type, assembler, note)
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
