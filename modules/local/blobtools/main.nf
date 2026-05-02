include { assemblerPublishDir } from '../../../lib/workflow_helpers.nf'

process BLOBTOOLS {
    tag "${sra}:${srr}:${assembler}"
    publishDir { assemblerPublishDir(sra, srr, assembler) }, mode: 'copy', overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(assembly_fasta), path(blast), path(assembly_bam), path(assembly_csi)
    path taxdump

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(assembly_fasta), path("blobtools.csv"), optional: true, emit: blobtable
    tuple val(sra), val(srr), val(read_type), val(assembler), path("blobtools*.svg"),                                                               optional: true, emit: blobplots
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("FAIL.note"),                         optional: true, emit: note

    script:
    """
    run_blobtools.sh \\
      --assembly "${assembly_fasta}" \\
      --bam "${assembly_bam}" \\
      --csi "${assembly_csi}" \\
      --blast "${blast}" \\
      --taxdump "${taxdump}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries}
    """
}
