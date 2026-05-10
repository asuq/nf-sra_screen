include { assemblerPublishDir } from '../../../lib/workflow_helpers.nf'

process BINETTE {
    tag "${sra}:${srr}:${assembler}"
    label 'binning'
    publishDir { "${assemblerPublishDir(sra, srr, assembler)}/binning" },
      mode: 'copy',
      overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler),
          path(assembly_fasta),
          path(contig2bin_maps)
    path checkm2_db

    output:
    tuple val(sra), val(srr), val(read_type), val(assembler), path("binette"),                                               emit: bins
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("binette.note"), emit: note

    script:
    def mapArg = contig2bin_maps.collect { contig2bin_map -> contig2bin_map.toString() }.join(',')
    """
    run_binette.sh \\
      --assembly "${assembly_fasta}" \\
      --bin-maps "${mapArg}" \\
      --checkm2-db "${checkm2_db}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries}
    """

    stub:
    """
    mkdir -p binette
    : > binette/final_contig_to_bin.tsv
    : > binette.note
    """
}
