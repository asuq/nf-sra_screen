include { assemblerPublishDir } from '../../../lib/workflow_helpers.nf'

process MAP_TO_ASSEMBLY {
    tag "${sra}:${srr}:${assembler}"
    publishDir { assemblerPublishDir(sra, srr, assembler) },
      mode: 'copy',
      overwrite: true,
      saveAs: { filename ->
        filename in [
          "assembly.bam.csi",
          "FAIL.note"
        ] ? filename : null
      }

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(assembly_fasta), path(reads)

    output:
    tuple val(sra), val(srr), val(read_type), val(assembler), path("assembly.bam"), path("assembly.bam.csi"),              optional: true, emit: mapped
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("FAIL.note"), optional: true, emit: note

    script:
    """
    run_map_to_assembly.sh \\
      --read-type "${read_type}" \\
      --assembly "${assembly_fasta}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries} \\
      --reads ${reads}
    """

    stub:
    """
    : > assembly.bam
    : > assembly.bam.csi
    """
}
