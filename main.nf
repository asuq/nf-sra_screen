#!/usr/bin/env nextflow
nextflow.enable.dsl=2


//-- Help Message ---------------------------------------------------------------

def helpMessage() {
  log.info """
  ===============================
  nf-sra_screen
  Nextflow pipeline for screening SRA genomes
  Version: 0.2.0
  Author : Akito Shima (ASUQ)
  Email: asuq.4096@gmail.com
  ===============================
  Usage: nextflow run main.nf [parameters]

  Required parameters:
    --sra           Path to sra.csv (header: sra)
    --taxa          Taxa for extraction (e.g., phylum, genus)
    --taxdump       Path to taxdump database folder
    --gtdb_ncbi_map Path to folder with GTDB-NCBI mapping Excel files
    --sandpiper_db  Path to Sandpiper database folder
    --singlem_db    Path to SingleM database folder
    --uniprot_db    Path to Uniprot database (.dmnd)

  Optional parameters:
    --help          Show this help message
    --outdir        Output directory (default: ./output)
    --max_retries   Maximum number of retries for each process (default: 3)
  """.stripIndent()
}


def missingParametersError() {
    log.error "Missing input parameters"
    helpMessage()
    error "Please provide all required parameters: --sra, --taxa, --taxdump, --gtdb_ncbi_map, --sandpiper_db, --singlem_db, and --uniprot_db"
}


//-- Processes -----------------------------------------------------------------

process VALIDATE_TAXA {
    input:
    path taxa_file
    path taxdump
    path gtdb_ncbi_map

    output:
    path("validated_taxa.csv"), emit: valid_taxa

    script:
    """
    jsonify_taxdump.py ${taxdump} \\
      && validate_taxa.py --taxa ${taxa_file} --taxdump ${taxdump} \\
        --gtdb-map --ncbi_to_gtdb ${gtdb_ncbi_map} --out validated_taxa.csv
    """
}


process DOWNLOAD_SRA_METADATA {
    tag "${sra}"
    publishDir "${params.outdir}/metadata/", mode: 'copy', overwrite: true

    input:
    val sra
    path valid_taxa

    output:
    tuple val(sra), path("${sra}.filtered.csv"), emit: filtered_sra
    tuple val(sra), path("${sra}.skipped.csv"),  emit: skipped_sra

    script:
    """
    # Retrieve metadata
    iseq -m -i "${sra}"

    # Extract SRR to screen
    filter_sra.sh "${sra}"
    """
}

process SANDPIPER {
    tag "${sra}:${srr}"
    publishDir "${params.outdir}/${sra}/${srr}/", mode: 'copy', overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler)
    path valid_taxa
    path sandpiper_db_ch

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("sandpiper_decision.txt"),   emit: decision
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("FAIL.note"), optional:true, emit: note
    tuple val(sra), val(srr), path("sandpiper_report.txt"),                                                optional:true, emit: sandpiper_report
    tuple val(sra), val(srr), path("sandpiper_output.tsv"),                                                optional:true, emit: sandpiper_summary

    script:
    """
    function hard_fail() {
      local msg="\$1"
      echo "\$msg" >&2
      if [[ ${task.attempt} -lt ${params.max_retries} ]]; then
        exit 1
      fi
      echo "\$msg" > FAIL.note
      echo "RUN_SINGLEM" > sandpiper_decision.txt
      exit 0
    }

    # Extract validated phyla
    awk -F',' '
      NR==1 {
        for (i = 1; i <= NF; i++) if (\$i == "gtdb_phylum") { c = i; break }
        next
      }
      c && \$c != "" { seen[\$c] = 1 }

      END {
        print "phyla"
        for (v in seen) print v
      }
      ' "${valid_taxa}" > phyla_to_check.txt \\
    || hard_fail "Sandpiper: failed to parse validated_taxa"

    # Lookup SRR in Sandpiper DB
    sandpiper_lookup.sh "${srr}" "${sandpiper_db_ch}" > sandpiper_report.txt \\
      || hard_fail "Sandpiper lookup failed"

    # Check if sandpiper result has taxa interested
    if check_singlem_phyla.py -i sandpiper_report.txt -p phyla_to_check.txt -o sandpiper_output.tsv; then
      # at least one target phylum found
      echo "PASS" > sandpiper_decision.txt
      exit 0
    else
      rc=\$?
      case "\$rc" in
        2)  # no target phyla
            echo "NEGATIVE" > sandpiper_decision.txt
            echo "No target phyla detected by Sandpiper" > FAIL.note
            exit 0;;
        1|*) # internal error – don’t kill the sample, just fallback
            hard_fail "Sandpiper phylum check error";;
      esac
    fi
    """
}


process DOWNLOAD_SRR {
    tag "${sra}:${srr}"

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), val(sandpiper)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), val(sandpiper), path("*.f*q*"), path("assembler.txt"), optional:true, emit: reads
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), val(sandpiper), path("FAIL.note"),                     optional:true, emit: note

    script:
    """
    # Download sequence data
    if ! iseq -Q -g -t ${task.cpus} -p 8 -i "${srr}"; then
      if [[ ${task.attempt} -lt ${params.max_retries} ]]; then
        exit 1
      fi
      echo "Download raw data failed" > FAIL.note
      rm -f *.f*q* "${srr}"  # remove sra and fastq files that didn't download properly
      exit 0
    fi

    # Pacbio assembler check
    final_asm="${assembler}"
    if [[ "${platform}" == "PACBIO_SMRT" && ( -z "${assembler}" || "${assembler}" == "unknown" ) ]]; then
      echo "Checking PacBio reads to determine assembler"
      if zcat -f *.f*q* 2>/dev/null \
        | awk 'NR%4==1{ h=tolower(\$0); if (h ~ /\\/ccs(\s|\$)/) { found=1; exit } } END{ exit(!found) }'
      then
        final_asm="long_hifi"
      else
        final_asm="long_pacbio"
      fi
      echo "\${final_asm}" > assembler.txt
    else
      touch assembler.txt
    fi
    """
}


process SINGLEM {
    tag "${sra}:${srr}"
    publishDir "${params.outdir}/${sra}/${srr}/", mode: 'copy', overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), val(sandpiper), path(reads)
    path valid_taxa
    path singlem_db_ch

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("reads_ok/*.f*q*"), optional:true, emit: reads
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("FAIL.note"),       optional:true, emit: note
    tuple val(sra), val(srr), path("singlem_taxonomic_profile*"),                                                optional:true, emit: singlem_summary
    tuple val(sra), val(srr), path("singlem_output.tsv"),                                                        optional:true, emit: singlem_phyla_check

    script:
    """
    function fail() {
      local msg="\$1"
      echo "\$msg" >&2
      if [[ ${task.attempt} -lt ${params.max_retries} ]]; then
        exit 1
      fi
      echo "\$msg" > FAIL.note
      exit 0
    }

    # Check if sandpiper already passed
    if [[ "${sandpiper}" == "PASS" ]]; then
      echo "Reads already passed Sandpiper check; skipping SingleM"
      mkdir -p reads_ok
      cp -v *.f*q* reads_ok/
      exit 0
    fi

    # Extract validated phyla
    awk -F',' '
      NR==1 {
        for (i = 1; i <= NF; i++) if (\$i == "gtdb_phylum") { c = i; break }
        next
      }
      c && \$c != "" { seen[\$c] = 1 }

      END {
        print "phyla"
        for (v in seen) print v
      }
      ' "${valid_taxa}" > phyla_to_check.txt \\
    || fail "SingleM: failed to parse validated_taxa"

    # Run singlem pipe
    R1=\$(ls *_1.fastq.gz *_R1*.fastq.gz 2>/dev/null | head -n1 || true)
    R2=\$(ls *_2.fastq.gz *_R2*.fastq.gz 2>/dev/null | head -n1 || true)
    if [[ -n "\$R1" && -n "\$R2" ]]; then
      singlem pipe -1 \$R1 -2 \$R2 \\
        --taxonomic-profile singlem_taxonomic_profile.tsv \\
        --taxonomic-profile-krona singlem_taxonomic_profile_krona \\
        --metapackage "${singlem_db_ch}" --threads ${task.cpus} \\
        || fail "SingleM pipe failed"
    else
      singlem pipe -1 ${reads} \\
        --taxonomic-profile singlem_taxonomic_profile.tsv \\
        --taxonomic-profile-krona singlem_taxonomic_profile_krona \\
        --metapackage "${singlem_db_ch}" --threads ${task.cpus} \\
        || fail "SingleM pipe failed"
    fi

    if [[ ! -s singlem_taxonomic_profile.tsv ]]; then
      fail "SingleM produced empty taxonomic profile"
    fi

    # Summarise singlem output
    singlem summarise \\
      --input-taxonomic-profile singlem_taxonomic_profile.tsv \\
      --output-species-by-site-relative-abundance-prefix singlem_taxonomic_profile_summary \\
      || fail "SingleM summarise failed"

    # Check for presence of target phyla
    check_singlem_phyla.sh singlem_taxonomic_profile.tsv phyla_to_check.txt singlem_output.tsv \\
      ${task.attempt} ${params.max_retries}
    """
}


process METASPADES {
    tag "${sra}:${srr}"
    publishDir "${params.outdir}/${sra}/${srr}/", mode: 'copy', overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path(reads)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("assembly.fasta"), optional:true, emit: assembly_fasta
    tuple val(sra), val(srr), path("assembly.gfa"),                                                             optional:true, emit: assembly_graph
    tuple val(sra), val(srr), path("spades.log"),                                                               optional:true, emit: assembly_log
    tuple val(sra), val(srr), path("assembly.bam"), path("assembly.bam.csi"),                                   optional:true, emit: assembly_bam
    tuple val(sra), val(srr), path("fastp.html"),                                                               optional:true, emit: fastp_html
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("FAIL.note"),      optional:true, emit: note

    script:
    """
    R1=\$(ls *_1.fastq.gz *_R1*.fastq.gz 2>/dev/null | head -n1 || true)
    R2=\$(ls *_2.fastq.gz *_R2*.fastq.gz 2>/dev/null | head -n1 || true)
    if [[ -z "\$R1" || -z "\$R2" ]]; then
      if [[ ${task.attempt} -lt ${params.max_retries} ]]; then exit 1; fi
      echo "Assembly failed: paired-end reads not found" > FAIL.note; exit 0
    fi

    # Run fastp
    if ! fastp --in1 "\${R1}" --in2 "\${R2}" --out1 "${srr}_fastp_R1.fastq.gz" --out2 "${srr}_fastp_R2.fastq.gz" \\
        --thread ${task.cpus} --length_required 50 --detect_adapter_for_pe --html fastp.html; then

      if [[ ${task.attempt} -lt ${params.max_retries} ]]; then exit 1; fi
      echo "Assembly failed at fastp" > FAIL.note; exit 0
    fi

    # Run SPAdes
    if ! spades.py -1 "${srr}_fastp_R1.fastq.gz" -2 "${srr}_fastp_R2.fastq.gz" -o '.' \\
      -k 21,33,55,77,99,119,127 --meta --threads ${task.cpus} --memory ${task.memory.toGiga()}; then

      if [[ ${task.attempt} -lt ${params.max_retries} ]]; then exit 1; fi
      echo "Assembly failed at SPAdes" > FAIL.note; exit 0
    fi

    # Rename outputs
    if [ -f scaffolds.fasta ]; then
      mv -v scaffolds.fasta assembly.fasta ;
    else
      mv -v contigs.fasta assembly.fasta ;
    fi
    mv -v assembly_graph_with_scaffolds.gfa assembly.gfa

    # Run bowtie2
    if ! ( bowtie2-build -f -q --threads ${task.cpus} assembly.fasta assembly_index \\
          && bowtie2 -q --reorder --threads ${task.cpus} --time --met-stderr --met 10 \\
            -x assembly_index -1 "${srr}_fastp_R1.fastq.gz" -2 "${srr}_fastp_R2.fastq.gz" \\
            | samtools sort --output-fmt BAM -@ ${task.cpus} -o assembly.bam \\
          && samtools index -c -o assembly.bam.csi -@ ${task.cpus} assembly.bam ); then

      if [[ ${task.attempt} -lt ${params.max_retries} ]]; then exit 1; fi
      echo "Assembly failed at bowtie2 mapping/indexing" > FAIL.note; exit 0
    fi
    """
}


process METAFLYE_NANO {
    tag "${sra}:${srr}"
    publishDir "${params.outdir}/${sra}/${srr}/", mode: 'copy', overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path(reads)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("assembly.fasta"), optional:true, emit: assembly_fasta
    tuple val(sra), val(srr), path("assembly.gfa"),                                                             optional:true, emit: assembly_graph
    tuple val(sra), val(srr), path("flye.log"),                                                                 optional:true, emit: assembly_log
    tuple val(sra), val(srr), path("assembly.bam"), path("assembly.bam.csi"),                                   optional:true, emit: assembly_bam
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("FAIL.note"),      optional:true, emit: note

    script:
    """
    # Run metaFlye
    if ! flye --nano-raw ${reads} --threads ${task.cpus} --scaffold --out-dir '.' --meta; then
      if [[ ${task.attempt} -lt ${params.max_retries} ]]; then exit 1; fi
      echo "Assembly failed at metaFlye (ONT)" > FAIL.note; exit 0
    fi

    # Run minimap2
    if ! ( minimap2 -ax map-ont -t ${task.cpus} assembly.fasta ${reads} \\
      | samtools sort --output-fmt BAM -@ ${task.cpus} -o assembly.bam \\
      && samtools index -c -o assembly.bam.csi -@ ${task.cpus} assembly.bam ); then

      if [[ ${task.attempt} -lt ${params.max_retries} ]]; then exit 1; fi
      echo "Assembly failed at mapping/indexing (ONT)" > FAIL.note; exit 0
    fi

    # Rename outputs
    mv -v assembly_graph.gfa assembly.gfa
    """
}


process METAFLYE_PACBIO {
    tag "${sra}:${srr}"
    publishDir "${params.outdir}/${sra}/${srr}/", mode: 'copy', overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path(reads)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("assembly.fasta"), optional:true, emit: assembly_fasta
    tuple val(sra), val(srr), path("assembly.gfa"),                                                             optional:true, emit: assembly_graph
    tuple val(sra), val(srr), path("flye.log"),                                                                 optional:true, emit: assembly_log
    tuple val(sra), val(srr), path("assembly.bam"), path("assembly.bam.csi"),                                   optional:true, emit: assembly_bam
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("FAIL.note"),      optional:true, emit: note

    script:
    """
    # Run metaFlye
    if ! flye --pacbio-raw ${reads} --threads ${task.cpus} --scaffold --out-dir '.' --meta; then
      if [[ ${task.attempt} -lt ${params.max_retries} ]]; then exit 1; fi
      echo "Assembly failed at metaFlye (PacBio)" > FAIL.note; exit 0
    fi

    # Run minimap2
    if ! ( minimap2 -ax map-pb -t ${task.cpus} assembly.fasta ${reads} \\
      | samtools sort --output-fmt BAM -@ ${task.cpus} -o assembly.bam \\
      && samtools index -c -o assembly.bam.csi -@ ${task.cpus} assembly.bam ); then

      if [[ ${task.attempt} -lt ${params.max_retries} ]]; then exit 1; fi
      echo "Assembly failed at mapping/indexing (PacBio)" > FAIL.note; exit 0
    fi

    # Rename outputs
    mv -v assembly_graph.gfa assembly.gfa
    """
}


process MYLOASM {
    tag "${sra}:${srr}"
    publishDir "${params.outdir}/${sra}/${srr}/", mode: 'copy', overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path(reads)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("assembly.fasta"), optional:true, emit: assembly_fasta
    tuple val(sra), val(srr), path("assembly.gfa"),                                                             optional:true, emit: assembly_graph
    tuple val(sra), val(srr), path("myloasm_*.log"),                                                            optional:true, emit: assembly_log
    tuple val(sra), val(srr), path("assembly.bam"), path("assembly.bam.csi"),                                   optional:true, emit: assembly_bam
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("FAIL.note"),      optional:true, emit: note

    script:
    """
    # Run myloasm
    if ! myloasm ${reads} -o . -t ${task.cpus} --hifi; then
      if [[ ${task.attempt} -lt ${params.max_retries} ]]; then exit 1; fi
      echo "Assembly failed at myloasm" > FAIL.note; exit 0
    fi

    # Rename outputs
    mv -v assembly_primary.fa assembly.fasta
    mv -v final_contig_graph.gfa assembly.gfa

    # Run minimap2
    if ! ( minimap2 -ax map-hifi -t ${task.cpus} assembly.fasta ${reads} \\
      | samtools sort --output-fmt BAM -@ ${task.cpus} -o assembly.bam \\
      && samtools index -c -o assembly.bam.csi -@ ${task.cpus} assembly.bam ); then

      if [[ ${task.attempt} -lt ${params.max_retries} ]]; then exit 1; fi
      echo "Assembly failed at mapping/indexing (HiFi)" > FAIL.note; exit 0
    fi
    """
}


process DIAMOND {
    tag "${sra}:${srr}"
    publishDir "${params.outdir}/${sra}/${srr}/", mode: 'copy', overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path(assembly_fasta)
    path uniprot_db

    output:
    tuple val(sra), val(srr), path("assembly_vs_uniprot.tsv"),                                             optional:true, emit: blast
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("FAIL.note"), optional:true, emit: note

    script:
    """
    if ! diamond blastx --sensitive --query "${assembly_fasta}" \\
        --out "assembly_vs_uniprot.tsv" --db "${uniprot_db}" \\
        --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \\
        --verbose --threads ${task.cpus} --evalue 1e-25 --max-target-seqs 5; then
      if [[ ${task.attempt} -lt ${params.max_retries} ]]; then exit 1; fi
      echo "Diamond failed" > FAIL.note; exit 0
    fi
    """
}


process BLOBTOOLS {
    tag "${sra}:${srr}"
    publishDir "${params.outdir}/${sra}/${srr}/", mode: 'copy', overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path(assembly_fasta), path(blast), path(assembly_bam), path(assembly_csi)
    path taxdump

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path(assembly_fasta), path("blobtools.csv"), optional:true, emit: blobtable
    tuple val(sra), val(srr), path("blobtools*.svg"),                                                                                optional:true, emit: blobplots
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("FAIL.note"),                           optional:true, emit: note

    script:
    """
    if ! blobtools create --fasta "${assembly_fasta}" --cov "${assembly_bam}" \\
      --hits "${blast}" --taxdump "${taxdump}" --threads  ${task.cpus} 'blobtools'; then
      if [[ ${task.attempt} -lt ${params.max_retries} ]]; then exit 1; fi
      echo "Blobtools failed at create" > FAIL.note; exit 0
    fi

    if ! blobtools filter --table 'blobtools.tsv' \\
      --table-fields gc,length,ncount,assembly_cov,bestsumorder_superkingdom,bestsumorder_kingdom,bestsumorder_phylum,bestsumorder_class,bestsumorder_order,bestsumorder_family,bestsumorder_genus,bestsumorder_species 'blobtools'; then
      if [[ ${task.attempt} -lt ${params.max_retries} ]]; then exit 1; fi
      echo "Blobtools failed at filter" > FAIL.note; exit 0
    fi

    blobtools view --format svg --plot 'blobtools' || {
      echo "blobtools view failed; leaving note and continuing." >&2
      echo "blobtools view failed" > blobtools/_view_failed.txt
    }

    # Convert tsv file to csv file
    header=\$(head -n 1 'blobtools.tsv' | cut -f2- | sed 's/bestsumorder_//g' | sed "s/assembly_cov/coverage/" | tr '\t' ',')
    data=\$(tail -n +2 'blobtools.tsv' | cut -f2- | tr '\t' ',')
    { echo "\${header}"; echo "\${data}"; } > 'blobtools.csv'
    """
}


process EXTRACT_TAXA {
    tag "${sra}:${srr}"
    publishDir "${params.outdir}/${sra}/${srr}/", mode: 'copy', overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path(assembly_fasta), path(blobtable)
    path(taxa)
    path(taxdump)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("summary.csv"), optional:true, emit: summary
    tuple val(sra), val(srr), path("*.ids.csv"),                                                             optional:true, emit: extracted_ids
    tuple val(sra), val(srr), path("*.fasta"),                                                               optional:true, emit: extracted_fasta
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("FAIL.note"),   optional:true, emit: note


    script:
    """
    if ! extract_records.py --blobtable "${blobtable}" \\
      --fasta "${assembly_fasta}" --taxa "${taxa}" --taxdump "${taxdump}"; then
      if [[ ${task.attempt} -lt ${params.max_retries} ]]; then exit 1; fi
      echo "Extraction failed" > FAIL.note; exit 0
    fi
    """
}


process LOG_FAILED_PROCESS {
  tag "${sra}"

  input:
  tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), val(note)

  output:
  tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("empty_summary.csv"), val(note), emit: skipped_rows

  script:
  """
  echo 'rank,taxa,n_identifiers,output_ids_csv,output_fasta' > empty_summary.csv
  """
}


process APPEND_SUMMARY {
    tag "${sra}:${srr}"

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path(summary_csv), val(note)
    val outdir

    output:
    path("summary.tsv"), optional: true, emit: global_summary

    script:
    """
    OUT_TSV="${outdir}/summary.tsv"
    mkdir -p "${outdir}"

    # counts = comma-joined n_identifiers in the order of rows in summary.csv
    COUNTS=\$(awk -F, 'NR>1{print \$3}' "${summary_csv}" | paste -sd, -)

    # note left blank for successful path
    LINE="${sra}\t${srr}\t${platform}\t${model}\t${strategy}\t${assembler}\t\${COUNTS}\t${note}"
    {
      flock 200
      if [[ ! -s "\$OUT_TSV" ]]; then
        echo -e "sra\tsrr\tplatform\tmodel\tstrategy\tassembler\tcounts\tnote" > "\$OUT_TSV"
      fi
      echo -e "\$LINE" >> "\$OUT_TSV"
    } 200> "\$OUT_TSV.lock"

    ln -sf "\$OUT_TSV" .
    """
}


//-- Workflow ------------------------------------------------------------------
workflow {
  if (params.help) {
    helpMessage()
    exit 0
  }

  if (!params.sra || !params.uniprot_db || !params.taxa \
      || !params.taxdump || !params.gtdb_ncbi_map \
      || !params.sandpiper_db || !params.singlem_db) {
    missingParametersError()
    exit 1
  }

  // Channel setup
  outdir = file(params.outdir).toAbsolutePath().toString()
  sra_ch = channel.fromPath(params.sra, checkIfExists: true)
                  .splitCsv(header: true, strip: true)
                  .map { row -> row.sra.trim() }
                  .filter { it }
                  .distinct()
  taxa_ch         = channel.value( file(params.taxa) )
  taxdump_ch      = channel.value( file(params.taxdump) )
  gtdb_ncbi_map   = channel.value( file(params.gtdb_ncbi_map) )
  singlem_db_ch   = channel.value( file(params.singlem_db) )
  sandpiper_db_ch = channel.value( file(params.sandpiper_db) )
  uniprot_db_ch   = channel.value( file(params.uniprot_db) )

  // Validate taxa
  validated_taxa = VALIDATE_TAXA(taxa_ch, taxdump_ch, gtdb_ncbi_map).valid_taxa

  // Step 1: extract metadata & filter SRR
  sra_metadata = DOWNLOAD_SRA_METADATA(sra_ch, validated_taxa)
  filtered_srr = sra_metadata.filtered_sra
  // skipped_srr.view { sra, csvfile -> "SKIPPED SRA: ${sra} (see ${csvfile})" }

  // Build nested channels per CSV, then flatten
  srr_ch = filtered_srr.map { sra, csvfile -> file(csvfile)}
            .splitCsv(header: true, strip: true)
            .map { row ->
                def sra = (row.accession ?: '').trim()
                def srr = (row.run_accession ?: '').trim()
                def platform = (row.instrument_platform ?: '').trim()
                def model = (row.instrument_model ?: '').trim()
                def strategy = (row.library_strategy ?: '').trim()
                def assembler = (row.assembler ?: '').trim()
                [sra, srr, platform, model, strategy, assembler]
            }
            .filter { it[1] }  // Ensure run_accession (SRR) is not empty
            .distinct()

  // srr_ch.view { acc, srr, platform, model, strategy, asm ->
  //       "DEBUG: ${acc}\t${srr}\t${platform}\t${model}\t${strategy}\t${asm}"
  // }

  // Step 2: check if srr has sandpiper results
  sandpiper = SANDPIPER(srr_ch, validated_taxa, sandpiper_db_ch)
  // decision: NEGATIVE / RUN_SINGLEM / PASS
  sandpiper_decision_ch = sandpiper.decision.map { sra, srr, platform, model, strategy, assembler, dec_path ->
      def decision = file(dec_path).text.trim()
      tuple(sra, srr, platform, model, strategy, assembler, decision)
  }

  srr_prescreened = sandpiper_decision_ch
    .filter { sra, srr, platform, model, strategy, assembler, decision ->
        decision == 'PASS' || decision == 'RUN_SINGLEM'
    }

  // Step 3: download SRR reads
  download_srr = DOWNLOAD_SRR(srr_prescreened)
  srr_reads = download_srr.reads.map { sra, srr, platform, model, strategy, asm, sandpiper_dec, reads, asm_txt ->
    def fixedAsm = file(asm_txt)?.text?.trim() ?: asm
    tuple(sra, srr, platform, model, strategy, fixedAsm, sandpiper_dec, reads)
  }

  // srr_reads.view { acc, srr, platform, model, strategy, asm, reads ->
  //   "${acc}\t${srr}\t${platform}\t${model}\t${strategy}\t${asm}\t${reads}"
  // }

  // Step 4: run SingleM to screen reads
  singlem = SINGLEM(srr_reads, validated_taxa, singlem_db_ch)
  singlem_reads = singlem.reads

  // Step 5: assemble reads
  short_ch       = singlem_reads.filter { sra, srr, platform, model, strategy, assembler, reads -> assembler.equalsIgnoreCase('short') }
  long_nano_ch   = singlem_reads.filter { sra, srr, platform, model, strategy, assembler, reads -> assembler.equalsIgnoreCase('long_nano') }
  long_pacbio_ch = singlem_reads.filter { sra, srr, platform, model, strategy, assembler, reads -> assembler.equalsIgnoreCase('long_pacbio') }
  long_hifi_ch   = singlem_reads.filter { sra, srr, platform, model, strategy, assembler, reads -> assembler.equalsIgnoreCase('long_hifi') }

  spades_asm     = METASPADES(short_ch)
  flyenano_asm   = METAFLYE_NANO(long_nano_ch)
  flyepacbio_asm = METAFLYE_PACBIO(long_pacbio_ch)
  hifimeta_asm   = MYLOASM(long_hifi_ch)

  // Step 6: run DIAMOND
  asm_fasta_ch = channel.empty()
                        .mix(spades_asm.assembly_fasta)
                        .mix(flyenano_asm.assembly_fasta)
                        .mix(flyepacbio_asm.assembly_fasta)
                        .mix(hifimeta_asm.assembly_fasta)


  diamond = DIAMOND(asm_fasta_ch, uniprot_db_ch)
  // diamond.blast.view { sra, srr, blast ->
  //   "DIAMOND: ${sra}\t${srr}\t${blast}"
  // }

  // Step 7. run BlobTools
  // Merge all BAM+Bai streams
  bam_ch = channel.empty()
                  .mix(spades_asm.assembly_bam)
                  .mix(flyenano_asm.assembly_bam)
                  .mix(flyepacbio_asm.assembly_bam)
                  .mix(hifimeta_asm.assembly_bam)

  // Key every stream by (sra,srr,assembler)
  fasta_by  = asm_fasta_ch.map   { sra, srr, platform, model, strategy, assembler, fasta -> tuple([sra,srr], [platform,model,strategy,assembler,fasta]) }
  blast_by  = diamond.blast.map  { sra, srr, hits  -> tuple([sra,srr], hits) }
  bam_by    = bam_ch.map         { sra, srr, bam, bai -> tuple([sra,srr], [bam,bai]) }

  // Join (fasta * diamond) then * bam
  fasta_blast = fasta_by.join(blast_by)
  fasta_blast_bam = fasta_blast.join(bam_by)

  // Unkey + call BlobTools
  blobtools_in = fasta_blast_bam.map { key, fasta, hits, pair ->
    def (sra, srr) = key
    def (platform, model, strategy, assembler, assembly) = fasta
    def (bam, bai) = pair
    tuple(sra, srr, platform, model, strategy, assembler, assembly, hits, bam, bai)
  }

  blobtools = BLOBTOOLS(blobtools_in, taxdump_ch)
  // blobtools.blobtable.view { sra, srr, platform, model, strategy, assembler, assembly, tbl -> "BLOBTOOLS:\t${sra}\t${srr}\t${platform}\t${model}\t${strategy}\t${assembler}\t${assembly}\t${tbl}" }

  // Step 8: extract taxa
  taxa_extraction = EXTRACT_TAXA(blobtools.blobtable, validated_taxa, taxdump_ch)

  // Step 9: append to global summary
  skipped_srr = sra_metadata.skipped_sra
    .map { sra, csvfile -> file(csvfile) }
    .splitCsv(header: true, strip: true)
    .map { row ->
      def sra   = (row.accession ?: '').trim()
      def srr   = (row.run_accession ?: '').trim()
      def plat  = (row.instrument_platform ?: '').trim()
      def model = (row.instrument_model ?: '').trim()
      def strat = (row.library_strategy ?: '').trim()
      def note  = "did not match the criteria: ${(row.skip_reason ?: '').trim()}"
      // assembler is empty for skipped rows
      tuple(sra, srr, plat, model, strat, '', note)
    }
    .filter { it[1] } // keep only rows with srr

  // Collect all errors
  errors = channel.empty()
                  .mix(sandpiper.note)
                  .mix(download_srr.note)
                  .mix(singlem.note)
                  .mix(spades_asm.note)
                  .mix(flyenano_asm.note)
                  .mix(flyepacbio_asm.note)
                  .mix(hifimeta_asm.note)
                  .mix(diamond.note)
                  .mix(blobtools.note)
                  .mix(taxa_extraction.note)
                  .map { sra, srr, platform, model, strategy, assembler, note_path ->
                    def note = file(note_path).text.trim()
                    tuple(sra, srr, platform, model, strategy, assembler, note)
                  }
                  .mix(skipped_srr)


  failed_sra = LOG_FAILED_PROCESS(errors)
  succeeded_sra = taxa_extraction.summary.map { sra, srr, platform, model, strategy, assembler, summary_csv ->
                                                tuple(sra, srr, platform, model, strategy, assembler, summary_csv, '')
                                              }

  // Combine succeeded and failed
  summary = channel.empty()
                   .mix(succeeded_sra)
                   .mix(failed_sra)

  APPEND_SUMMARY(summary, outdir)
}
