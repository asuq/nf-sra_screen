process CREATE_EMPTY_SUMMARY {
  tag "${sra}:${srr}:${read_type}:${assembler}"
  label 'run_local'

  input:
  tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), val(note)

  output:
  tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("empty_summary.csv"), val(note), emit: skipped_rows

  script:
  """
  echo 'rank,ncbi_taxa,n_contigs,output_ids_csv,output_fasta' > empty_summary.csv
  """
}
