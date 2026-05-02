process RESOLVE_SRR_METADATA {
    tag "${sra}:${srr}"

    input:
    tuple val(sra), val(srr), path(assembly_fasta)

    output:
    tuple val(sra), val(srr), path("resolved.tsv"), path(assembly_fasta), optional: true, emit: resolved
    tuple val(sra), val(srr), val('UNKNOWN'), val('UNKNOWN'), val('UNKNOWN'), val('UNKNOWN'), val('provided'), path("FAIL.note"), optional: true, emit: note

    script:
    def resolveScript = file("${workflow.projectDir}/bin/run_resolve_srr_metadata.sh").toAbsolutePath()
    """
    ${resolveScript} \\
      --sample "${sra}" \\
      --srr "${srr}" \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries}
    """

    stub:
    """
    printf 'ILLUMINA\tNovaSeq 6000\tWGS\tshort\n' > resolved.tsv
    """
}

