include { assemblerPublishDir } from '../../../lib/workflow_helpers.nf'

process SEMIBIN {
    tag "${sra}:${srr}:${assembler}"
    label 'binning'
    publishDir { "${assemblerPublishDir(sra, srr, assembler)}/binning" },
      mode: 'copy',
      overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(assembly_fasta), path(assembly_bam), path(assembly_csi)
    path uniprot_db

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler),
          val("semibin"), path("semibin"), path("semibin.contig2bin.tsv"), path("semibin.note"),             emit: result

    script:
    """
    run_semibin.sh \\
      --assembly "${assembly_fasta}" \\
      --bam "${assembly_bam}" \\
      --diamond-db "${uniprot_db}" \\
      --read-type "${read_type}" \\
      --environment "${params.semibin_environment}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries}
    """
}
