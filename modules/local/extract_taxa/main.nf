include { assemblerPublishDir } from '../../../lib/workflow_helpers.nf'

process EXTRACT_TAXA {
    tag "${sra}:${srr}:${assembler}"
    publishDir { assemblerPublishDir(sra, srr, assembler) }, mode: 'copy', overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(assembly_fasta), path(blobtable)
    path valid_taxa
    path taxdump

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("summary.csv"), optional: true, emit: summary
    tuple val(sra), val(srr), val(read_type), val(assembler), path("*.ids.csv"),                                            optional: true, emit: extracted_ids
    tuple val(sra), val(srr), val(read_type), val(assembler), path("*.fasta"),                                              optional: true, emit: extracted_fasta
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("FAIL.note"),  optional: true, emit: note

    script:
    """
    run_extract_taxa.sh \\
      --blobtable "${blobtable}" \\
      --fasta "${assembly_fasta}" \\
      --taxa "${valid_taxa}" \\
      --taxdump "${taxdump}" \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries}
    """
}
