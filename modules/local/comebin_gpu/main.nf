include { assemblerPublishDir } from '../../../lib/workflow_helpers.nf'

process COMEBIN_GPU {
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
          val("comebin"), path("comebin"), path("comebin.contig2bin.tsv"), path("FAIL.note"),                emit: result

    script:
    """
    run_comebin_nf.sh \\
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
    : > FAIL.note
    """
}
