include { assemblerPublishDir } from '../../../lib/workflow_helpers.nf'

process FASTP_SHORT {
    tag "${sra}:${srr}:${assembler}"
    publishDir { assemblerPublishDir(sra, srr, assembler) },
      mode: 'copy',
      overwrite: true,
      saveAs: { filename ->
        filename in [
          "fastp.html",
          "fastp.json",
          "FAIL.note"
        ] ? filename : null
      }

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(reads)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("*_fastp_R*.fastq.gz"), optional: true, emit: reads
    tuple val(sra), val(srr), val(read_type), val(assembler), path("fastp.html"),                                                        optional: true, emit: html
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("FAIL.note"),               optional: true, emit: note

    script:
    """
    run_fastp_short.sh \\
      --srr "${srr}" \\
      --assembler "${assembler}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries} \\
      --reads ${reads}
    """
}
