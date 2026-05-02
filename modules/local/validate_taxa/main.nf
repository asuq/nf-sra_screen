process VALIDATE_TAXA {
    input:
    path taxa_file
    path taxdump
    path gtdb_ncbi_map

    output:
    path("validated_taxa.csv"), emit: valid_taxa

    script:
    """
    jsonify_taxdump.py "${taxdump}" \\
    && validate_taxa.py --taxa "${taxa_file}" --taxdump "${taxdump}" \\
      --gtdb-map --ncbi_to_gtdb "${gtdb_ncbi_map}" --out validated_taxa.csv
    """
}
