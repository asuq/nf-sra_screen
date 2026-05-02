include { assemblerPublishDir } from '../../../lib/workflow_helpers.nf'

process UNICYCLER {
    tag "${sra}:${srr}:${assembler}"
    label 'assembly'
    publishDir { assemblerPublishDir(sra, srr, assembler) },
      mode: 'copy',
      overwrite: true,
      saveAs: { filename ->
        filename in [
          "assembly.fasta",
          "assembly.gfa",
          "spades.log",
          "FAIL.note"
        ] ? filename : null
      }

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(reads)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("assembly.fasta"), optional: true, emit: assembly_fasta
    tuple val(sra), val(srr), val(read_type), val(assembler), path("assembly.gfa"),                                            optional: true, emit: assembly_graph
    tuple val(sra), val(srr), val(read_type), val(assembler), path("spades.log"),                                              optional: true, emit: assembly_log
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("FAIL.note"),     optional: true, emit: note

    script:
    """
    # --reads should be the last argument and unquoted to capture all read files
    run_unicycler.sh \\
      --srr "${srr}" \\
      --cpus ${task.cpus} \\
      --memory-gb ${task.memory.toGiga()} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries} \\
      --reads ${reads}
    """
}
