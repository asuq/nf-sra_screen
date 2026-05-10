include { assemblerPublishDir } from '../../../lib/workflow_helpers.nf'

process LORBIN_GPU {
    tag "${sra}:${srr}:${assembler}"
    label 'gpu'
    label 'binning'
    publishDir { "${assemblerPublishDir(sra, srr, assembler)}/binning" },
      mode: 'copy',
      overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(assembly_fasta), path(assembly_bam), path(assembly_csi)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler),
          val("lorbin"), path("lorbin"), path("lorbin.contig2bin.tsv"), path("lorbin.note"),                  emit: result

    script:
    """
    run_lorbin.sh \\
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
