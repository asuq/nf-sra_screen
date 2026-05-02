process SINGLEM {
    tag "${sra}:${srr}"
    publishDir { "${params.outdir}/${sra}/${srr}/" },
      mode: 'copy',
      overwrite: true,
      saveAs: { filename ->
        if (filename == "singlem_output.tsv") return filename
        if (filename.startsWith("singlem_taxonomic_profile")) return filename
        if (filename == "FAIL.note") return filename
        return null
      }

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(sandpiper), path(reads)
    path valid_taxa
    path singlem_db

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), path("reads_ok/*.f*q*"),          optional: true, emit: reads
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(''), path("FAIL.note"),       optional: true, emit: note
    tuple val(sra), val(srr), path("singlem_taxonomic_profile*"),                                                optional: true, emit: singlem_summary
    tuple val(sra), val(srr), path("singlem_output.tsv"),                                                        optional: true, emit: singlem_phyla_check

    script:
    """
    # --reads should be the last argument and unquoted to capture all read files
    run_singlem.sh \\
      --sandpiper-decision "${sandpiper}" \\
      --valid-taxa "${valid_taxa}" \\
      --singlem-db "${singlem_db}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries} \\
      --reads ${reads}
    """
}
