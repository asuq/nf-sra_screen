include { assemblerPublishDir } from '../../../lib/workflow_helpers.nf'

process VAMB_GPU {
    tag "${sra}:${srr}:${assembler}"
    label 'gpu'
    publishDir { "${assemblerPublishDir(sra, srr, assembler)}/binning" },
      mode: 'copy',
      overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(assembly_fasta), path(assembly_bam), path(assembly_csi)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler),
          val("vamb"), path("vamb"), path("vamb.contig2bin.tsv"), path("vamb.note"),                         emit: result

    script:
    """
    run_vamb.sh \\
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
