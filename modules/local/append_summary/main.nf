process APPEND_SUMMARY {
    tag "${sra}:${srr}:${read_type}:${assembler}"

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(summary_csv), val(note)
    val outdir

    output:
    path("summary.tsv"), optional: true, emit: global_summary

    script:
    """
    run_append_summary.sh \\
      --outdir "${outdir}" \\
      --sra "${sra}" \\
      --srr "${srr}" \\
      --platform "${platform}" \\
      --model "${model}" \\
      --strategy "${strategy}" \\
      --read-type "${read_type}" \\
      --assembler "${assembler}" \\
      --summary-csv "${summary_csv}" \\
      --note "${note}"
    """
}
