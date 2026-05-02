include { assemblerPublishDir } from '../../../lib/workflow_helpers.nf'

process DIAMOND {
    tag "${sra}:${srr}:${assembler}"
    publishDir { assemblerPublishDir(sra, srr, assembler) }, mode: 'copy', overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(assembly_fasta)
    path uniprot_db

    output:
    tuple val(sra), val(srr), val(read_type), val(assembler), path("assembly_vs_uniprot.tsv"),                                     optional: true, emit: blast
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("FAIL.note"), optional: true, emit: note

    script:
    """
    run_diamond.sh \\
      --assembly "${assembly_fasta}" \\
      --db "${uniprot_db}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries}
    """
}
