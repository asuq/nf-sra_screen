include { assemblerPublishDir } from '../../../lib/workflow_helpers.nf'

process METAFLYE_HIFI {
    tag "${sra}:${srr}:${assembler}"
    label 'assembly'
    publishDir { assemblerPublishDir(sra, srr, assembler) },
      mode: 'copy',
      overwrite: true,
      saveAs: { filename ->
        filename in [
          "assembly.fasta",
          "assembly.gfa",
          "flye.log",
          "FAIL.note"
        ] ? filename : null
      }

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(reads)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("assembly.fasta"), optional: true, emit: assembly_fasta
    tuple val(sra), val(srr), val(read_type), val(assembler), path("assembly.gfa"),                                            optional: true, emit: assembly_graph
    tuple val(sra), val(srr), val(read_type), val(assembler), path("flye.log"),                                                optional: true, emit: assembly_log
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("FAIL.note"),     optional: true, emit: note

    script:
    """
    run_metaflye_hifi.sh \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries} \\
      --reads ${reads}
    """
}
