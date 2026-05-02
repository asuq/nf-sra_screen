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


include {
  validateBinningOptions
  binnersForSample
  binnerCsvContains
  binnerCsvSize
  normaliseBinningRow
  parseResolvedMetadata
} from './lib/workflow_helpers.nf'


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
    publishDir { "${params.outdir}/${sra}/${srr}/" },
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
    publishDir { "${params.outdir}/${sra}/${srr}/binning" },
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
    publishDir { "${params.outdir}/${sra}/${srr}/binning" },
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
    publishDir { "${params.outdir}/${sra}/${srr}/binning" },
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
    publishDir { "${params.outdir}/${sra}/${srr}/binning" },
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
    publishDir { "${params.outdir}/${sra}/${srr}/binning" },
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
    publishDir { "${params.outdir}/${sra}/${srr}/binning" },
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
    publishDir { "${params.outdir}/${sra}/${srr}/binning" },
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
    publishDir { "${params.outdir}/${sra}/${srr}/binning" },
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
    publishDir { "${params.outdir}/${sra}/${srr}/binning" },
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
    publishDir { "${params.outdir}/${sra}/${srr}/binning" },
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
    def mapArg = contig2bin_maps.collect { contig2bin_map -> contig2bin_map.toString() }.join(',')
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
    publishDir { "${params.outdir}/${sra}/${srr}/binning" },
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
    def mapArg = contig2bin_maps.collect { contig2bin_map -> contig2bin_map.toString() }.join(',')
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
    main:
    if (params.help) {
      helpMessage()
      exit 0
    }

    if (!params.binning_tsv || !params.uniprot_db) {
      missingParametersError()
    }

    def binning_options = validateBinningOptions()
    def selected_binners = binning_options.binners
    def selected_refiners = binning_options.refiners
    def use_gpu_binners = binning_options.gpu

    if (use_gpu_binners) {
      def cpu_only_binners = selected_binners.findAll { binner -> binner in ['metabat', 'rosella', 'semibin'] }
      if (cpu_only_binners) {
        log.warn "GPU mode requested; CPU-only binner(s) will remain on CPU: ${cpu_only_binners.join(', ')}"
      }
    }

    def outdir = file(params.outdir ?: './output').toAbsolutePath().toString()
    def uniprot_db_ch = channel.value(file(params.uniprot_db, checkIfExists: true))

    def binning_rows_channel = channel.fromPath(params.binning_tsv, checkIfExists: true)
      .splitCsv(header: true, sep: '\t', strip: true)
      .map { row -> normaliseBinningRow(row) }

    def local_binning_rows_channel = binning_rows_channel
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

    def srr_binning_rows_channel = binning_rows_channel
      .filter { row -> row.mode == 'srr' }
      .map { row ->
        tuple(row.sra, row.srr, row.assembly_fasta)
      }

    def resolved_srr_out = RESOLVE_SRR_METADATA(srr_binning_rows_channel)
    def resolved_srr_rows_channel = resolved_srr_out.resolved.map { sra, srr, resolved_tsv, assembly_fasta ->
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

    def downloaded_srr_out = DOWNLOAD_SRR(resolved_srr_rows_channel)
    def downloaded_binning_rows_channel = downloaded_srr_out.reads.map { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, reads, asm_txt ->
      def fixedReadType = file(asm_txt).text.trim() ?: read_type
      tuple(sra, srr, platform, model, strategy, fixedReadType, assembler, assembly_fasta, reads)
    }

    def mapping_input_rows_channel = channel.empty()
      .mix(local_binning_rows_channel)
      .mix(downloaded_binning_rows_channel)

    def short_mapping_rows_channel = mapping_input_rows_channel
      .filter { _sra, _srr, _platform, _model, _strategy, read_type, _assembler, _assembly_fasta, _reads ->
        read_type.equalsIgnoreCase('short')
      }
    def nanopore_mapping_rows_channel = mapping_input_rows_channel
      .filter { _sra, _srr, _platform, _model, _strategy, read_type, _assembler, _assembly_fasta, _reads ->
        read_type.equalsIgnoreCase('nanopore')
      }
    def pacbio_mapping_rows_channel = mapping_input_rows_channel
      .filter { _sra, _srr, _platform, _model, _strategy, read_type, _assembler, _assembly_fasta, _reads ->
        read_type.equalsIgnoreCase('pacbio')
      }
    def hifi_mapping_rows_channel = mapping_input_rows_channel
      .filter { _sra, _srr, _platform, _model, _strategy, read_type, _assembler, _assembly_fasta, _reads ->
        read_type.equalsIgnoreCase('hifi')
      }

    def short_mapping_out = MAP_SHORT(short_mapping_rows_channel)
    def nanopore_mapping_out = MAP_NANO(nanopore_mapping_rows_channel)
    def pacbio_mapping_out = MAP_PACBIO(pacbio_mapping_rows_channel)
    def hifi_mapping_out = MAP_HIFI(hifi_mapping_rows_channel)

    def mapped_bam_channel = channel.empty()
      .mix(short_mapping_out.mapped)
      .mix(nanopore_mapping_out.mapped)
      .mix(pacbio_mapping_out.mapped)
      .mix(hifi_mapping_out.mapped)

    def assembly_fasta_by_sample_key = mapping_input_rows_channel.map { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, _reads ->
      tuple([sra, srr, read_type, assembler], [platform, model, strategy, assembly_fasta])
    }

    def mapped_bam_by_sample_key = mapped_bam_channel.map { sra, srr, _platform, _model, _strategy, read_type, assembler, assembly_bam, assembly_csi ->
      tuple([sra, srr, read_type, assembler], [assembly_bam, assembly_csi])
    }

    def mapped_assembly_rows_channel = assembly_fasta_by_sample_key
      .join(mapped_bam_by_sample_key)
      .map { sample_join_key, assembly_payload, bam_payload ->
        def (sra, srr, read_type, assembler) = sample_join_key
        def (platform, model, strategy, assembly_fasta) = assembly_payload
        def (assembly_bam, assembly_csi) = bam_payload
        tuple(sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi)
      }

    def mapping_failure_note_channel = channel.empty()
      .mix(resolved_srr_out.note)
      .mix(downloaded_srr_out.note)
      .mix(short_mapping_out.note)
      .mix(nanopore_mapping_out.note)
      .mix(pacbio_mapping_out.note)
      .mix(hifi_mapping_out.note)

    def planned_binning_input = mapped_assembly_rows_channel
      .map { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi ->
        def sample_binner_csv = binnersForSample(binning_options, sra, srr, read_type).join(',')
        tuple(sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi, sample_binner_csv)
      }

    def binner_plan_by_sample = planned_binning_input
      .map { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi, sample_binner_csv ->
        tuple([sra, srr, read_type, assembler], binnerCsvSize(sample_binner_csv))
      }

    def metabat_input = planned_binning_input
      .filter { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi, sample_binner_csv ->
        binnerCsvContains(sample_binner_csv, 'metabat')
      }
      .map { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi, sample_binner_csv ->
        tuple(sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi)
      }

    def comebin_input = planned_binning_input
      .filter { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi, sample_binner_csv ->
        binnerCsvContains(sample_binner_csv, 'comebin')
      }
      .map { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi, sample_binner_csv ->
        tuple(sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi)
      }

    def vamb_input = planned_binning_input
      .filter { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi, sample_binner_csv ->
        binnerCsvContains(sample_binner_csv, 'vamb')
      }
      .map { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi, sample_binner_csv ->
        tuple(sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi)
      }

    def lorbin_input = planned_binning_input
      .filter { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi, sample_binner_csv ->
        binnerCsvContains(sample_binner_csv, 'lorbin')
      }
      .map { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi, sample_binner_csv ->
        tuple(sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi)
      }

    def semibin_input = planned_binning_input
      .filter { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi, sample_binner_csv ->
        binnerCsvContains(sample_binner_csv, 'semibin')
      }
      .map { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi, sample_binner_csv ->
        tuple(sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi)
      }

    def rosella_input = planned_binning_input
      .filter { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi, sample_binner_csv ->
        binnerCsvContains(sample_binner_csv, 'rosella')
      }
      .map { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi, sample_binner_csv ->
        tuple(sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, assembly_bam, assembly_csi)
      }

    def metabat_results = channel.empty()
    def comebin_results = channel.empty()
    def vamb_results = channel.empty()
    def lorbin_results = channel.empty()
    def semibin_results = channel.empty()
    def rosella_results = channel.empty()

    if ('metabat' in selected_binners) {
      metabat_results = METABAT(metabat_input).result
    }
    if ('comebin' in selected_binners) {
      if (use_gpu_binners) {
        comebin_results = COMEBIN_GPU(comebin_input).result
      }
      else {
        comebin_results = COMEBIN(comebin_input).result
      }
    }
    if ('vamb' in selected_binners) {
      if (use_gpu_binners) {
        vamb_results = VAMB_GPU(vamb_input).result
      }
      else {
        vamb_results = VAMB(vamb_input).result
      }
    }
    if ('lorbin' in selected_binners) {
      if (use_gpu_binners) {
        lorbin_results = LORBIN_GPU(lorbin_input).result
      }
      else {
        lorbin_results = LORBIN(lorbin_input).result
      }
    }
    if ('semibin' in selected_binners) {
      semibin_results = SEMIBIN(semibin_input, uniprot_db_ch).result
    }
    if ('rosella' in selected_binners) {
      rosella_results = ROSELLA(rosella_input).result
    }

    def binner_results = channel.empty()
      .mix(metabat_results)
      .mix(comebin_results)
      .mix(vamb_results)
      .mix(lorbin_results)
      .mix(semibin_results)
      .mix(rosella_results)

    def dastool_base_by = mapped_assembly_rows_channel.map { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, _bam, _csi ->
      tuple([sra, srr, read_type, assembler], [platform, model, strategy, assembly_fasta])
    }

    def binner_maps_by_sample = binner_results
      .map { sra, srr, platform, model, strategy, read_type, assembler, tool, bin_dir, contig2bin, note_path ->
        tuple([sra, srr, read_type, assembler], [tool, contig2bin])
      }
      .combine(binner_plan_by_sample, by: 0)
      .map { sample_join_key, binner_map_payload, expectedCount -> tuple(groupKey(sample_join_key, expectedCount as int), binner_map_payload) }
      .groupTuple()
      .map { grouped_sample_key, binner_map_entries -> tuple(grouped_sample_key.target, binner_map_entries) }

    def dastool_join = dastool_base_by.join(binner_maps_by_sample)

    def dastool_in = dastool_join.map { sample_join_key, assembly_payload, binner_map_entries ->
      def (sra, srr, read_type, assembler) = sample_join_key
      def (platform, model, strategy, assembly_fasta) = assembly_payload
      def contig2bin_maps = binner_map_entries.collect { binner_map_entry -> binner_map_entry[1] }
      tuple(
        sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta,
        contig2bin_maps
      )
    }

    def dastool_note_entries = channel.empty()
    if ('dastool' in selected_refiners) {
      dastool_note_entries = DASTOOL(dastool_in).note.map { sra, srr, platform, model, strategy, read_type, assembler, note_path ->
        tuple(sra, srr, platform, model, strategy, read_type, assembler, 'dastool', note_path)
      }
    }

    def binette_note_entries = channel.empty()
    if ('binette' in selected_refiners) {
      def checkm2_db_ch = channel.value(file(params.checkm2_db, checkIfExists: true))
      binette_note_entries = BINETTE(dastool_in, checkm2_db_ch).note.map { sra, srr, platform, model, strategy, read_type, assembler, note_path ->
        tuple(sra, srr, platform, model, strategy, read_type, assembler, 'binette', note_path)
      }
    }

    def successful_mapping_rows_channel = mapped_assembly_rows_channel
      .map { sra, srr, platform, model, strategy, read_type, assembler, _assembly_fasta, _bam, _csi ->
        tuple(sra, srr, platform, model, strategy, read_type, assembler, '')
      }
      .distinct()

    def empty_summary_out = CREATE_EMPTY_SUMMARY(successful_mapping_rows_channel).skipped_rows
    def succeeded_rows = empty_summary_out.map { sra, srr, platform, model, strategy, read_type, assembler, summary_csv, _note ->
      tuple(sra, srr, platform, model, strategy, read_type, assembler, summary_csv, '')
    }

    def binner_note_entries = binner_results.map { sra, srr, platform, model, strategy, read_type, assembler, tool, bin_dir, contig2bin, note_path ->
      tuple(sra, srr, platform, model, strategy, read_type, assembler, tool, note_path)
    }

    def binning_notes_by_sample_channel = channel.empty()
      .mix(binner_note_entries)
      .mix(dastool_note_entries)
      .mix(binette_note_entries)
      .map { sra, srr, platform, model, strategy, read_type, assembler, tool, note_path ->
        tuple([sra, srr, platform, model, strategy, read_type, assembler], [tool, note_path])
      }
      .groupTuple()
      .map { sample_key, note_entries ->
        def (sra, srr, platform, model, strategy, read_type, assembler) = sample_key
        def note_text = note_entries
          .sort { left, right -> left[0] <=> right[0] }
          .collect { note_entry ->
            def note_path = note_entry[1]
            def text = file(note_path as String).text.trim()
            text ? "${note_entry[0]}: ${text}" : null
          }
          .findAll { token -> token }
          .join('; ')
        tuple(sra, srr, platform, model, strategy, read_type, assembler, note_text)
      }

    def final_success_rows_channel = succeeded_rows
      .map { sra, srr, platform, model, strategy, read_type, assembler, summary_csv, note ->
        tuple([sra, srr, platform, model, strategy, read_type, assembler], [summary_csv, note])
      }
      .join(binning_notes_by_sample_channel.map { sra, srr, platform, model, strategy, read_type, assembler, note ->
        tuple([sra, srr, platform, model, strategy, read_type, assembler], note)
      })
      .map { sample_key, summary_payload, binning_note ->
        def (sra, srr, platform, model, strategy, read_type, assembler) = sample_key
        def summary_csv = summary_payload[0]
        def base_note = summary_payload[1] ?: ''
        def extra_note = (binning_note ?: '').trim()
        def final_note = base_note
        if (extra_note) {
          final_note = final_note ? "${final_note}; ${extra_note}" : extra_note
        }
        tuple(sra, srr, platform, model, strategy, read_type, assembler, summary_csv, final_note)
      }

    def mapping_errors = mapping_failure_note_channel.map { sra, srr, platform, model, strategy, read_type, assembler, note_path ->
      def note = file(note_path).text.trim()
      tuple(sra, srr, platform, model, strategy, read_type, assembler, note)
    }

    def failed_rows = CREATE_EMPTY_FAILURE_SUMMARY(mapping_errors)

    def summary_rows_channel = channel.empty()
      .mix(final_success_rows_channel)
      .mix(failed_rows)

    APPEND_SUMMARY(summary_rows_channel, outdir)

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
