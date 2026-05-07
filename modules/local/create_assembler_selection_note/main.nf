process CREATE_ASSEMBLER_SELECTION_NOTE {
  tag "${sra}:${srr}:${read_type}"
  label 'run_local'

  input:
  tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), val(note)

  output:
  tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("FAIL.note"), emit: note

  script:
  """
  printf '%s\n' "${note}" > FAIL.note
  """
}
