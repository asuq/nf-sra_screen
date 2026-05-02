include { assemblerPublishDir } from '../../../lib/workflow_helpers.nf'

process ROSELLA {
    tag "${sra}:${srr}:${assembler}"
    label 'binning'
    publishDir { "${assemblerPublishDir(sra, srr, assembler)}/binning" },
      mode: 'copy',
      overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(assembly_fasta), path(assembly_bam), path(assembly_csi)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler),
          val("rosella"), path("rosella"), path("rosella.contig2bin.tsv"), path("rosella.note"),             emit: result

    script:
    """
    run_rosella.sh \\
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
