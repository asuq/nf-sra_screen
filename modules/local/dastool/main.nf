include { assemblerPublishDir } from '../../../lib/workflow_helpers.nf'

process DASTOOL {
    tag "${sra}:${srr}:${assembler}"
    label 'binning'
    publishDir { "${assemblerPublishDir(sra, srr, assembler)}/binning" },
      mode: 'copy',
      overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler),
          path(assembly_fasta),
          path(contig2bin_maps)

    output:
    tuple val(sra), val(srr), val(read_type), val(assembler), path("dastool"),                                               emit: bins
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("dastool.note"), emit: note

    script:
    def mapArg = contig2bin_maps.collect { contig2bin_map -> contig2bin_map.toString() }.join(',')
    """
    run_dastool.sh \\
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
