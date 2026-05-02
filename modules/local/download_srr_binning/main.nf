process DOWNLOAD_SRR {
    tag "${sra}:${srr}"
    publishDir { "${params.outdir}/${sra}/${srr}/" },
      mode: 'copy',
      overwrite: true,
      saveAs: { filename ->
        filename in ["FAIL.note"] ? filename : null
      }

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(assembly_fasta)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(assembly_fasta), path("*.f*q*"), path("assembler.txt"), optional: true, emit: reads
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("FAIL.note"), optional: true, emit: note

    script:
    """
    run_download_srr.sh \\
      --srr "${srr}" \\
      --platform "${platform}" \\
      --read-type "${read_type}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries}
    """

    stub:
    """
    : > "${srr}_1.fastq.gz"
    printf '%s\n' "${read_type}" > assembler.txt
    """
}

