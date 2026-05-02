process DOWNLOAD_SRA_METADATA {
    tag "${sra}"
    publishDir { "${params.outdir}/metadata/${sra}/" },
        mode: 'copy',
        overwrite: true,
        saveAs: { filename ->
            if (filename == "FAIL.note") {
                return "${sra}.FAIL.note"
            }
            return filename
        }

    input:
    val sra

    output:
    tuple val(sra), path("${sra}.filtered.csv"),                                             optional: true, emit: filtered_sra
    tuple val(sra), path("${sra}.skipped.csv"),                                              optional: true, emit: skipped_sra
    tuple val(sra), val(''), val(''), val(''), val(''), val(''), val(''), path("FAIL.note"), optional: true, emit: note

    script:
    """
    run_download_metadata.sh \\
      --sra "${sra}" \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries}
    """
}
