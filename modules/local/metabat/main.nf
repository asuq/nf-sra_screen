include { assemblerPublishDir } from '../../../lib/workflow_helpers.nf'

process METABAT {
    tag "${sra}:${srr}:${assembler}"
    label 'binning'
    publishDir { "${assemblerPublishDir(sra, srr, assembler)}/binning" },
      mode: 'copy',
      overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(assembly_fasta), path(assembly_bam), path(assembly_csi)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler),
          val("metabat"), path("metabat"), path("metabat.contig2bin.tsv"), path("metabat.note"),             emit: result

    script:
    """
    run_metabat.sh \\
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
