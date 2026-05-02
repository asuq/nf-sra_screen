process SANDPIPER {
    tag "${sra}:${srr}"
    publishDir { "${params.outdir}/${sra}/${srr}/" }, mode: 'copy', overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type)
    path valid_taxa
    path sandpiper_db

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), path("sandpiper_decision.txt"),            emit: decision
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(''), path("FAIL.note"), optional: true, emit: note
    tuple val(sra), val(srr), path("sandpiper_report.txt"),                                                optional: true, emit: sandpiper_report
    tuple val(sra), val(srr), path("sandpiper_output.tsv"),                                                optional: true, emit: sandpiper_summary

    script:
    """
    run_sandpiper.sh \\
      --srr "${srr}" \\
      --valid-taxa "${valid_taxa}" \\
      --db "${sandpiper_db}" \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries}
    """
}
