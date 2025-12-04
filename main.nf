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
    --taxa          Path to taxa.csv for extraction (header: rank, taxa)
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
    jsonify_taxdump.py "${taxdump}" \\
    && validate_taxa.py --taxa "${taxa_file}" --taxdump "${taxdump}" \\
      --gtdb-map --ncbi_to_gtdb "${gtdb_ncbi_map}" --out validated_taxa.csv
    """
}


process DOWNLOAD_SRA_METADATA {
    tag "${sra}"
    publishDir "${params.outdir}/metadata/${sra}/",
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
    path valid_taxa

    output:
    tuple val(sra), path("${sra}.filtered.csv"),                                    optional: true, emit: filtered_sra
    tuple val(sra), path("${sra}.skipped.csv"),                                     optional: true, emit: skipped_sra
    tuple val(sra), val(''), val(''), val(''), val(''), val(''), path("FAIL.note"), optional: true, emit: note

    script:
    """
    run_download_metadata.sh \\
      --sra "${sra}" \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries}
    """
}


process SANDPIPER {
    tag "${sra}:${srr}"
    publishDir "${params.outdir}/${sra}/${srr}/", mode: 'copy', overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler)
    path valid_taxa
    path sandpiper_db

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("sandpiper_decision.txt"),    emit: decision
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("FAIL.note"), optional: true, emit: note
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


process DOWNLOAD_SRR {
    tag "${sra}:${srr}"
    publishDir "${params.outdir}/${sra}/${srr}/",
      mode: 'copy',
      overwrite: true,
      saveAs: { filename ->
        filename in ["FAIL.note"] ? filename : null
      }

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), val(sandpiper)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), val(sandpiper), path("*.f*q*"), path("assembler.txt"), optional: true, emit: reads
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("FAIL.note"),                                     optional: true, emit: note

    script:
    """
    run_download_srr.sh \\
      --srr "${srr}" \\
      --platform "${platform}" \\
      --assembler "${assembler}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries}
    """
}


process SINGLEM {
    tag "${sra}:${srr}"
    publishDir "${params.outdir}/${sra}/${srr}/",
      mode: 'copy',
      overwrite: true,
      saveAs: { filename ->
        if (filename == "singlem_output.tsv") return filename
        if (filename.startsWith("singlem_taxonomic_profile")) return filename
        if (filename == "FAIL.note") return filename
        return null
      }

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), val(sandpiper), path(reads)
    path valid_taxa
    path singlem_db

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("reads_ok/*.f*q*"), optional: true, emit: reads
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("FAIL.note"),       optional: true, emit: note
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


process METASPADES {
    tag "${sra}:${srr}"
    label 'assembly'
    publishDir "${params.outdir}/${sra}/${srr}/",
      mode: 'copy',
      overwrite: true,
      saveAs: { filename ->
        filename in [
          "assembly.fasta",
          "assembly.gfa",
          "spades.log",
          "assembly.bam.csi",
          "fastp.html",
          "FAIL.note"
        ] ? filename : null
      }

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path(reads)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("assembly.fasta"), optional: true, emit: assembly_fasta
    tuple val(sra), val(srr), path("assembly.gfa"),                                                             optional: true, emit: assembly_graph
    tuple val(sra), val(srr), path("spades.log"),                                                               optional: true, emit: assembly_log
    tuple val(sra), val(srr), path("assembly.bam"), path("assembly.bam.csi"),                                   optional: true, emit: assembly_bam
    tuple val(sra), val(srr), path("fastp.html"),                                                               optional: true, emit: fastp_html
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("FAIL.note"),      optional: true, emit: note

    script:
    """
    # --reads should be the last argument and unquoted to capture all read files
    run_metaspades.sh \\
      --srr "${srr}" \\
      --cpus ${task.cpus} \\
      --memory-gb ${task.memory.toGiga()} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries} \\
      --reads ${reads}
    """
}


process METAFLYE_NANO {
    tag "${sra}:${srr}"
    label 'assembly'
    publishDir "${params.outdir}/${sra}/${srr}/",
      mode: 'copy',
      overwrite: true,
      saveAs: { filename ->
        filename in [
          "assembly.fasta",
          "assembly.gfa",
          "flye.log",
          "assembly.bam.csi",
          "FAIL.note"
        ] ? filename : null
      }

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path(reads)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("assembly.fasta"), optional: true, emit: assembly_fasta
    tuple val(sra), val(srr), path("assembly.gfa"),                                                             optional: true, emit: assembly_graph
    tuple val(sra), val(srr), path("flye.log"),                                                                 optional: true, emit: assembly_log
    tuple val(sra), val(srr), path("assembly.bam"), path("assembly.bam.csi"),                                   optional: true, emit: assembly_bam
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("FAIL.note"),      optional: true, emit: note

    script:
    """
    run_metaflye_nano.sh \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries} \\
      --reads ${reads}
    """
}


process METAFLYE_PACBIO {
    tag "${sra}:${srr}"
    label 'assembly'
    publishDir "${params.outdir}/${sra}/${srr}/",
      mode: 'copy',
      overwrite: true,
      saveAs: { filename ->
        filename in [
          "assembly.fasta",
          "assembly.gfa",
          "flye.log",
          "assembly.bam.csi",
          "FAIL.note"
        ] ? filename : null
      }

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path(reads)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("assembly.fasta"), optional: true, emit: assembly_fasta
    tuple val(sra), val(srr), path("assembly.gfa"),                                                             optional: true, emit: assembly_graph
    tuple val(sra), val(srr), path("flye.log"),                                                                 optional: true, emit: assembly_log
    tuple val(sra), val(srr), path("assembly.bam"), path("assembly.bam.csi"),                                   optional: true, emit: assembly_bam
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("FAIL.note"),      optional: true, emit: note

    script:
    """
    run_metaflye_pacbio.sh \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries} \\
      --reads ${reads}
    """
}


process MYLOASM {
    tag "${sra}:${srr}"
    label 'assembly'
    publishDir "${params.outdir}/${sra}/${srr}/",
      mode: 'copy',
      overwrite: true,
      saveAs: { filename ->
        if (filename == "assembly.fasta") return filename
        if (filename == "assembly.gfa") return filename
        if (filename == "assembly.bam.csi") return filename
        if (filename.startsWith("myloasm_") && filename.endsWith(".log")) return filename
        if (filename == "FAIL.note") return filename
        return null
      }

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path(reads)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("assembly.fasta"), optional: true, emit: assembly_fasta
    tuple val(sra), val(srr), path("assembly.gfa"),                                                             optional: true, emit: assembly_graph
    tuple val(sra), val(srr), path("myloasm_*.log"),                                                            optional: true, emit: assembly_log
    tuple val(sra), val(srr), path("assembly.bam"), path("assembly.bam.csi"),                                   optional: true, emit: assembly_bam
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("FAIL.note"),      optional: true, emit: note

    script:
    """
    run_myloasm_hifi.sh \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries} \\
      --reads ${reads}
    """
}


process DIAMOND {
    tag "${sra}:${srr}"
    publishDir "${params.outdir}/${sra}/${srr}/", mode: 'copy', overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path(assembly_fasta)
    path uniprot_db

    output:
    tuple val(sra), val(srr), path("assembly_vs_uniprot.tsv"),                                             optional: true, emit: blast
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("FAIL.note"), optional: true, emit: note

    script:
    """
    run_diamond.sh \\
      --assembly "${assembly_fasta}" \\
      --db "${uniprot_db}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries}
    """
}


process BLOBTOOLS {
    tag "${sra}:${srr}"
    publishDir "${params.outdir}/${sra}/${srr}/", mode: 'copy', overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path(assembly_fasta), path(blast), path(assembly_bam), path(assembly_csi)
    path taxdump

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path(assembly_fasta), path("blobtools.csv"), optional: true, emit: blobtable
    tuple val(sra), val(srr), path("blobtools*.svg"),                                                                                optional: true, emit: blobplots
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("FAIL.note"),                           optional: true, emit: note

    script:
    """
    run_blobtools.sh \\
      --assembly "${assembly_fasta}" \\
      --bam "${assembly_bam}" \\
      --csi "${assembly_csi}" \\
      --blast "${blast}" \\
      --taxdump "${taxdump}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries}
    """
}


process EXTRACT_TAXA {
    tag "${sra}:${srr}"
    publishDir "${params.outdir}/${sra}/${srr}/", mode: 'copy', overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path(assembly_fasta), path(blobtable)
    path valid_taxa
    path taxdump

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("summary.csv"), optional: true, emit: summary
    tuple val(sra), val(srr), path("*.ids.csv"),                                                             optional: true, emit: extracted_ids
    tuple val(sra), val(srr), path("*.fasta"),                                                               optional: true, emit: extracted_fasta
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("FAIL.note"),   optional: true, emit: note

    script:
    """
    run_extract_taxa.sh \\
      --blobtable "${blobtable}" \\
      --fasta "${assembly_fasta}" \\
      --taxa "${valid_taxa}" \\
      --taxdump "${taxdump}" \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries}
    """
}


process METABAT {
    tag "${sra}:${srr}"
    label 'binning'
    publishDir "${params.outdir}/${sra}/${srr}/binning",
      mode: 'copy',
      overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path(assembly_fasta), path(assembly_bam), path(assembly_csi)

    output:
    tuple val(sra), val(srr), path("metabat"),                                                                emit: bins
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("metabat.note"), emit: note

    script:
    """
    run_metabat.sh \\
      --assembly "${assembly_fasta}" \\
      --bam "${assembly_bam}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries}
    """
}


process CONCOCT {
    tag "${sra}:${srr}"
    label 'binning'
    publishDir "${params.outdir}/${sra}/${srr}/binning",
      mode: 'copy',
      overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path(assembly_fasta), path(assembly_bam), path(assembly_csi)

    output:
    tuple val(sra), val(srr), path("concoct"),                                                                emit: bins
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("concoct.note"), emit: note

    script:
    """
    run_concoct.sh \\
      --assembly "${assembly_fasta}" \\
      --bam "${assembly_bam}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries}
    """
}


process SEMIBIN {
    tag "${sra}:${srr}"
    label 'binning'
    publishDir "${params.outdir}/${sra}/${srr}/binning",
      mode: 'copy',
      overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path(assembly_fasta), path(assembly_bam), path(assembly_csi)
    path uniprot_db

    output:
    tuple val(sra), val(srr), path("semibin"), path("semibin/contig_bins.tsv"),                               emit: bins
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("semibin.note"), emit: note

    script:
    """
    run_semibin.sh \\
      --assembly "${assembly_fasta}" \\
      --bam "${assembly_bam}" \\
      --diamond-db "${uniprot_db}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries}
    """
}


process ROSELLA {
    tag "${sra}:${srr}"
    label 'binning'
    publishDir "${params.outdir}/${sra}/${srr}/binning",
      mode: 'copy',
      overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path(assembly_fasta), path(assembly_bam), path(assembly_csi)

    output:
    tuple val(sra), val(srr), path("rosella"),                                                                emit: bins
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("rosella.note"), emit: note

    script:
    """
    run_rosella.sh \\
      --assembly "${assembly_fasta}" \\
      --bam "${assembly_bam}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries}
    """
}


process DASTOOL {
    tag "${sra}:${srr}"
    publishDir "${params.outdir}/${sra}/${srr}/binning",
      mode: 'copy',
      overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler),
          path(assembly_fasta),
          path(metabat_dir),
          path(concoct_dir),
          path(semibin_dir),
          path(semibin_contig_bins),
          path(rosella_dir)

    output:
    tuple val(sra), val(srr), path("dastool"),                                                                emit: bins
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("dastool.note"), emit: note

    script:
    def metaDir     = metabat_dir          ?: ''
    def concoctDir  = concoct_dir          ?: ''
    def semibinDir  = semibin_dir          ?: ''
    def semibinMap  = semibin_contig_bins  ?: ''
    def rosellaDir  = rosella_dir          ?: ''
    """
    run_dastool.sh \\
      --assembly "${assembly_fasta}" \\
      --metabat-dir "${metaDir}" \\
      --concoct-dir "${concoctDir}" \\
      --semibin-dir "${semibinDir}" \\
      --semibin-map "${semibinMap}" \\
      --rosella-dir "${rosellaDir}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries}
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
  echo 'rank,ncbi_taxa,n_contigs,output_ids_csv,output_fasta' > empty_summary.csv
  """
}


process BINNING_ERROR_SUMMARY {
    tag "${sra}:${srr}"

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler),
      path(metabat_note),
      path(concoct_note),
      path(semibin_note),
      path(rosella_note),
      path(dastool_note)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("binning_note.txt"), emit: note

    script:
    """
    : > binning_note.txt

    add_note() {
      local file="\$1"

      if [ -f "\$file" ]; then
        # Flatten file content to a single line
        note=\$(paste -sd' ' "\$file")

        # Skip if note is only whitespace
        if ! printf '%s' "\$note" | grep -q '[^[:space:]]'; then
          return 0
        fi

        # If we already wrote something, prepend a separator
        if [ -s binning_note.txt ]; then
          printf ';' >> binning_note.txt
        fi

        # Append the note as-is (no extra label; note already contains tool info)
        printf '%s' "\$note" >> binning_note.txt
      fi
    }

    add_note "${metabat_note}"
    add_note "${concoct_note}"
    add_note "${semibin_note}"
    add_note "${rosella_note}"
    add_note "${dastool_note}"

    # End the file with a single newline
    echo >> binning_note.txt
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
    run_append_summary.sh \\
      --outdir "${outdir}" \\
      --sra "${sra}" \\
      --srr "${srr}" \\
      --platform "${platform}" \\
      --model "${model}" \\
      --strategy "${strategy}" \\
      --assembler "${assembler}" \\
      --summary-csv "${summary_csv}" \\
      --note "${note}"
    """
}


//-- Workflow ------------------------------------------------------------------
workflow PRE_SCREENING {
  take:
    sra_ch
    validated_taxa
    sandpiper_db_ch
    singlem_db_ch

  main:
    // Step 1: extract metadata & filter SRR
    sra_metadata = DOWNLOAD_SRA_METADATA(sra_ch, validated_taxa)
    filtered_srr = sra_metadata.filtered_sra

    // Build nested channels per CSV, then flatten
    srr_ch = filtered_srr.map { sra, csvfile -> file(csvfile) }
              .splitCsv(header: true, strip: true)
              .map { row ->
                  def sra       = (row.accession ?: '').trim()
                  def srr       = (row.run_accession ?: '').trim()
                  def platform  = (row.instrument_platform ?: '').trim()
                  def model     = (row.instrument_model ?: '').trim()
                  def strategy  = (row.library_strategy ?: '').trim()
                  def assembler = (row.assembler ?: '').trim()
                  [sra, srr, platform, model, strategy, assembler]
              }
              .filter { it[1] }  // ensure SRR not empty
              .distinct()

    // Step 2: SANDPIPER prescreening
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

    // Step 4: SINGLEM prescreening
    singlem = SINGLEM(srr_reads, validated_taxa, singlem_db_ch)
    singlem_reads = singlem.reads

  emit:
    // For assembly
    singlem_reads         = singlem_reads

    // For summary
    sra_metadata_skipped  = sra_metadata.skipped_sra
    sra_metadata_note     = sra_metadata.note
    sandpiper_note        = sandpiper.note
    download_srr_note     = download_srr.note
    singlem_note          = singlem.note
}


workflow ASSEMBLY {
  take:
    singlem_reads
    validated_taxa
    uniprot_db_ch
    taxdump_ch

  main:
    // Step 5: assemble reads
    short_ch       = singlem_reads.filter { sra, srr, platform, model, strategy, assembler, reads -> assembler.equalsIgnoreCase('short') }
    long_nano_ch   = singlem_reads.filter { sra, srr, platform, model, strategy, assembler, reads -> assembler.equalsIgnoreCase('long_nano') }
    long_pacbio_ch = singlem_reads.filter { sra, srr, platform, model, strategy, assembler, reads -> assembler.equalsIgnoreCase('long_pacbio') }
    long_hifi_ch   = singlem_reads.filter { sra, srr, platform, model, strategy, assembler, reads -> assembler.equalsIgnoreCase('long_hifi') }

    spades_asm     = METASPADES(short_ch)
    flyenano_asm   = METAFLYE_NANO(long_nano_ch)
    flyepacbio_asm = METAFLYE_PACBIO(long_pacbio_ch)
    hifimeta_asm   = MYLOASM(long_hifi_ch)

    // Step 6: DIAMOND
    asm_fasta_ch = channel.empty()
                        .mix(spades_asm.assembly_fasta)
                        .mix(flyenano_asm.assembly_fasta)
                        .mix(flyepacbio_asm.assembly_fasta)
                        .mix(hifimeta_asm.assembly_fasta)

    diamond = DIAMOND(asm_fasta_ch, uniprot_db_ch)

    // BAMs for BlobTools and binning
    bam_src = channel.empty()
                  .mix(spades_asm.assembly_bam)
                  .mix(flyenano_asm.assembly_bam)
                  .mix(flyepacbio_asm.assembly_bam)
                  .mix(hifimeta_asm.assembly_bam)

    // Step 7: BlobTools
    // Key every stream by (sra,srr)
    fasta_by = asm_fasta_ch.map  { sra, srr, platform, model, strategy, assembler, fasta -> tuple([sra,srr], [platform,model,strategy,assembler,fasta]) }
    blast_by = diamond.blast.map { sra, srr, blast -> tuple([sra,srr], blast) }
    bam_by   = bam_src.map       { sra, srr, bam, csi -> tuple([sra,srr], [bam,csi]) }

    // Join (fasta * diamond) then * bam
    fasta_blast     = fasta_by.join(blast_by)
    fasta_blast_bam = fasta_blast.join(bam_by)

    // Unkey + call BlobTools
    blobtools_in = fasta_blast_bam.map { key, fasta, blast, pair ->
      def (sra, srr) = key
      def (platform, model, strategy, assembler, assembly) = fasta
      def (bam, csi) = pair
      tuple(sra, srr, platform, model, strategy, assembler, assembly, blast, bam, csi)
    }

    blobtools = BLOBTOOLS(blobtools_in, taxdump_ch)

    // Step 8: EXTRACT_TAXA
    taxa_extraction = EXTRACT_TAXA(blobtools.blobtable, validated_taxa, taxdump_ch)

    assembly_notes_ch = channel.empty()
                          .mix(spades_asm.note)
                          .mix(flyenano_asm.note)
                          .mix(flyepacbio_asm.note)
                          .mix(hifimeta_asm.note)

  emit:
    // For binning
    assembly_bam_all = bam_src
    blobtable        = blobtools.blobtable

    // For summary
    taxa_summary     = taxa_extraction.summary
    taxa_note        = taxa_extraction.note
    assembly_notes   = assembly_notes_ch
    diamond_note     = diamond.note
    blobtools_note   = blobtools.note
}


workflow BINNING {
  take:
    blobtable_ch
    assembly_bam_ch
    uniprot_db_ch

  main:
    // Build binning_input from blobtable + BAM
    // blobtable_ch: (sra, srr, platform, model, strategy, assembler, assembly_fasta, blobtable)
    // assembly_bam_ch: (sra, srr, bam, csi)
    bam_by = assembly_bam_ch.map { sra, srr, bam, csi -> tuple([sra, srr], [bam, csi]) }
    binning_base_by = blobtable_ch.map { sra, srr, platform, model, strategy, assembler, assembly_fasta, blobtable ->
      tuple([sra, srr], [platform, model, strategy, assembler, assembly_fasta])
    }
    binning_join = binning_base_by.join(bam_by)

    // binning_input: (sra, srr, platform, model, strategy, assembler, assembly_fasta, bam, csi)
    binning_input = binning_join.map { key, meta, bam_idx ->
      def (sra, srr) = key
      def (platform, model, strategy, assembler, assembly_fasta) = meta
      def (bam, csi) = bam_idx
      tuple(sra, srr, platform, model, strategy, assembler, assembly_fasta, bam, csi)
    }

    // Run individual binners
    metabat_binning = METABAT(binning_input)
    concoct_binning = CONCOCT(binning_input)
    semibin_binning = SEMIBIN(binning_input, uniprot_db_ch)
    rosella_binning = ROSELLA(binning_input)

    // Prepare DASTool input
    dastool_base_by = binning_input.map { sra, srr, platform, model, strategy, assembler, assembly_fasta, bam, csi ->
      tuple([sra,srr], [platform, model, strategy, assembler, assembly_fasta])
    }

    // Tag each entry with a label
    metabat_entries = metabat_binning.bins.map { sra, srr, metabat_dir ->
      tuple([sra, srr], metabat_dir)
    }
    concoct_entries = concoct_binning.bins.map { sra, srr, concoct_dir ->
      tuple([sra, srr], concoct_dir)
    }
    semibin_entries = semibin_binning.bins.map { sra, srr, semibin_dir, semibin_bins ->
      tuple([sra, srr], [semibin_dir, semibin_bins])
    }
    rosella_entries = rosella_binning.bins.map { sra, srr, rosella_dir ->
      tuple([sra, srr], rosella_dir)
    }

    // Group all tool-entries by sample
    dastool_join = dastool_base_by
      .join(metabat_entries)
      .join(concoct_entries)
      .join(semibin_entries)
      .join(rosella_entries)

    // Build DASTool input
    dastool_in = dastool_join.map { key, meta, metabat_dir, concoct_dir, semibin_pair, rosella_dir ->
      def (sra, srr) = key
      def (platform, model, strategy, assembler, assembly_fasta) = meta
      def (semibin_dir, semibin_bins) = semibin_pair
      tuple(
        sra, srr, platform, model, strategy, assembler, assembly_fasta,
        metabat_dir, concoct_dir, semibin_dir, semibin_bins, rosella_dir
      )
    }
    .filter { it != null }

    dastool_binning = DASTOOL(dastool_in)

    // Binning error aggregation
    def N_BIN_STATUS = 5

    metabat_status = metabat_binning.note.map { sra, srr, platform, model, strategy, assembler, note ->
      def key = groupKey([sra: sra, srr: srr, platform: platform, model: model, strategy: strategy, assembler: assembler], N_BIN_STATUS)
      tuple(key, ['metabat', note])
    }

    concoct_status = concoct_binning.note.map { sra, srr, platform, model, strategy, assembler, note ->
      def key = groupKey([sra: sra, srr: srr, platform: platform, model: model, strategy: strategy, assembler: assembler], N_BIN_STATUS)
      tuple(key, ['concoct', note])
    }

    semibin_status = semibin_binning.note.map { sra, srr, platform, model, strategy, assembler, note ->
      def key = groupKey([sra: sra, srr: srr, platform: platform, model: model, strategy: strategy, assembler: assembler], N_BIN_STATUS)
      tuple(key, ['semibin', note])
    }

    rosella_status = rosella_binning.note.map { sra, srr, platform, model, strategy, assembler, note ->
      def key = groupKey([sra: sra, srr: srr, platform: platform, model: model, strategy: strategy, assembler: assembler], N_BIN_STATUS)
      tuple(key, ['rosella', note])
    }

    dastool_status = dastool_binning.note.map { sra, srr, platform, model, strategy, assembler, note ->
      def key = groupKey([sra: sra, srr: srr, platform: platform, model: model, strategy: strategy, assembler: assembler], N_BIN_STATUS)
      tuple(key, ['dastool', note])
    }

    binning_mix = channel.empty()
      .mix(metabat_status)
      .mix(concoct_status)
      .mix(semibin_status)
      .mix(rosella_status)
      .mix(dastool_status)
      .groupTuple()

    binning_status_grouped = binning_mix
      .map { key, items ->
        def m = (Map) key.target
        def (sra, srr, platform, model, strategy, assembler) = [m.sra, m.srr, m.platform, m.model, m.strategy, m.assembler]
        def mp = items.collectEntries { tool, note -> [(tool): note] }
        tuple(sra, srr, platform, model, strategy, assembler, mp.metabat, mp.concoct, mp.semibin, mp.rosella, mp.dastool)
      }

    binning_error_summary = BINNING_ERROR_SUMMARY(binning_status_grouped)


  emit:
    binning_note = binning_error_summary.note
}


workflow SUMMARY {
  take:
    // from PRE_SCREENING
    sra_metadata_skipped
    sra_metadata_note
    sandpiper_note
    download_srr_note
    singlem_note

    // from ASSEMBLY
    assembly_notes
    diamond_note
    blobtools_note
    taxa_note
    taxa_summary

    // from BINNING
    binning_note

    // constant
    outdir

  main:
    // Turn skipped_sra CSVs into rows with a textual note
    skipped_srr = sra_metadata_skipped
      .map { sra, csvfile -> file(csvfile) }
      .splitCsv(header: true, strip: true)
      .map { row ->
        def sra       = (row.accession           ?: '').trim()
        def srr       = (row.run_accession       ?: '').trim()
        def platform  = (row.instrument_platform ?: '').trim()
        def model     = (row.instrument_model    ?: '').trim()
        def strategy  = (row.library_strategy    ?: '').trim()
        def note      = "did not match the criteria: ${(row.skip_reason ?: '').trim()}"
        // assembler is empty for skipped rows
        tuple(sra, srr, platform, model, strategy, '', note)
      }
      .filter { it[1] } // keep only rows with srr

    // Collect all non-binning errors
    errors = channel.empty()
                    .mix(sra_metadata_note)
                    .mix(sandpiper_note)
                    .mix(download_srr_note)
                    .mix(singlem_note)
                    .mix(assembly_notes)
                    .mix(diamond_note)
                    .mix(blobtools_note)
                    .mix(taxa_note)
                    .map { sra, srr, platform, model, strategy, assembler, note_path ->
                      def note = file(note_path).text.trim()
                      tuple(sra, srr, platform, model, strategy, assembler, note)
                    }
                    .mix(skipped_srr)

    failed_sra = LOG_FAILED_PROCESS(errors)

    // Successful samples: have a summary.csv (note starts empty)
    succeeded_sra = taxa_summary.map { sra, srr, platform, model, strategy, assembler, summary_csv ->
      tuple(sra, srr, platform, model, strategy, assembler, summary_csv, '')
    }

    // Combine succeeded_sra and binning_note by key
    def N_SUMMARY_ITEMS = 2
    succ_keyed = succeeded_sra.map { sra, srr, platform, model, strategy, assembler, summary_csv, base_note ->
      def keyMap = [sra: sra, srr: srr, platform: platform, model: model, strategy: strategy, assembler: assembler]
      def key = groupKey(keyMap, N_SUMMARY_ITEMS)
      tuple(key, [summary_csv, base_note])
    }

    binning_keyed = binning_note.map { sra, srr, platform, model, strategy, assembler, binning_note ->
      def keyMap = [sra: sra, srr: srr, platform: platform,
                    model: model, strategy: strategy, assembler: assembler]
      def key = groupKey(keyMap, N_SUMMARY_ITEMS)
      tuple(key, binning_note)
    }

    succ_and_bin = channel.empty()
                        .mix(succ_keyed)
                        .mix(binning_keyed)
                        .groupTuple()

    // Build final succeeded_sra with binning annotations in the note field
    succeeded_with_binning = succ_and_bin
      .map { key, values ->
        def m = (Map) key.target
        def (sra, srr, platform, model, strategy, assembler) =
          [m.sra, m.srr, m.platform, m.model, m.strategy, m.assembler]

        // Find the summary entry (List [summary_csv, base_note])
        def summaryEntry = values.find { it instanceof List && it.size() == 2 }
        if (!summaryEntry) {
          return null   // no summary => not a "successful" sample; ignore here
        }

        def summary_csv = summaryEntry[0]
        def base_note   = summaryEntry[1] ?: ''

        // All other values are binning notes (Strings)
        def binning_note_file = values.find { !(it instanceof List) }
        def binning_note_test = ''
        if (binning_note_file) {
          binning_note_test = file(binning_note_file as String).text.trim()
        }

        def final_note = base_note
        if (binning_note_test) {
          final_note = final_note ? "${final_note}; ${binning_note_test}" : binning_note_test
        }

        tuple(sra, srr, platform, model, strategy, assembler, summary_csv, final_note)
      }
      .filter { it != null }

    // Combine succeeded and failed
    summary = channel.empty()
                    .mix(succeeded_with_binning)
                    .mix(failed_sra)

    summary_result = APPEND_SUMMARY(summary, outdir)

  emit:
    global_summary = summary_result.global_summary
}


workflow {
    if (params.help) {
      helpMessage()
      exit 0
    }

    if (!params.sra || !params.uniprot_db || !params.taxa ||
        !params.taxdump || !params.gtdb_ncbi_map ||
        !params.sandpiper_db || !params.singlem_db) {
      missingParametersError()
    }

    // Channel setup
    outdir = file(params.outdir).toAbsolutePath().toString()
    sra_ch = channel.fromPath(params.sra, checkIfExists: true)
                    .splitCsv(header: true, strip: true)
                    .map { row -> row.sra.trim() }
                    .filter { it }
                    .distinct()
    taxa_ch           = channel.value( file(params.taxa) )
    taxdump_ch        = channel.value( file(params.taxdump) )
    gtdb_ncbi_map_ch  = channel.value( file(params.gtdb_ncbi_map) )
    singlem_db_ch     = channel.value( file(params.singlem_db) )
    sandpiper_db_ch   = channel.value( file(params.sandpiper_db) )
    uniprot_db_ch     = channel.value( file(params.uniprot_db) )

    // Step 0: validate taxa
    validated_taxa = VALIDATE_TAXA(taxa_ch, taxdump_ch, gtdb_ncbi_map_ch).valid_taxa

    // Step 1: PRE_SCREENING: DOWNLOAD_SRA_METADATA -> SANDPIPER -> DOWNLOAD_SRR -> SINGLEM
    pre = PRE_SCREENING(sra_ch, validated_taxa, sandpiper_db_ch, singlem_db_ch)

    // Step 2: ASSEMBLY: assemblers -> DIAMOND -> BLOBTOOLS -> EXTRACT_TAXA
    asm = ASSEMBLY(pre.singlem_reads, validated_taxa, uniprot_db_ch, taxdump_ch)

    // Step 3: BINNING: METABAT, CONCOCT, SEMIBIN, ROSELLA -> DASTOOL
    binning = BINNING(asm.blobtable, asm.assembly_bam_all, uniprot_db_ch)

    // Step 4: Build summary.tsv (including binning annotations)
    SUMMARY(
      // PRE_SCREENING: summary-related outputs
      pre.sra_metadata_skipped,
      pre.sra_metadata_note,
      pre.sandpiper_note,
      pre.download_srr_note,
      pre.singlem_note,

      // ASSEMBLY: summary-related outputs
      asm.assembly_notes,
      asm.diamond_note,
      asm.blobtools_note,
      asm.taxa_note,
      asm.taxa_summary,

      // BINNING: notes from each binner
      binning.binning_note,

      // constant
      outdir
    )
}

workflow.onComplete {
    def outdirPath  = file(params.outdir ?: './output').toAbsolutePath()
    def summaryFile = outdirPath.resolve('summary.tsv')
    def traceFile   = file("${workflow.launchDir}/execution-reports/trace.tsv").toAbsolutePath()
    def scriptFile  = file("${workflow.projectDir}/bin/annotate_summary_from_trace.py").toAbsolutePath()

    log.info "onComplete: summary.tsv -> ${summaryFile}"
    log.info "onComplete: trace.tsv   -> ${traceFile}"
    log.info "onComplete: annotator   -> ${scriptFile}"

    // Sanity checks
    if( !summaryFile.exists() ) {
      log.warn "onComplete: ${summaryFile} not found; skipping scheduler annotation"
      return
    }
    if( !traceFile.exists() ) {
      log.warn "onComplete: ${traceFile} not found; skipping scheduler annotation"
      return
    }

    // python3 annotate_summary_from_trace.py summary.tsv trace.tsv
    def cmd = [
      'python3',
      scriptFile.toString(),
      summaryFile.toString(),
      traceFile.toString()
    ]

    log.info "onComplete: running ${cmd.join(' ')}"

    // Run the script in the launch directory (where execution-report lives)
    def proc = new ProcessBuilder(cmd)
      .directory( workflow.launchDir.toFile() )
      .redirectError( java.lang.ProcessBuilder.Redirect.INHERIT )
      .redirectOutput( java.lang.ProcessBuilder.Redirect.INHERIT )
      .start()

    int rc = proc.waitFor()
    if( rc != 0 ) {
      log.warn "onComplete: annotator script exited with code ${rc}"
    }
    else {
      log.info "onComplete: summary.tsv successfully annotated with scheduler error information"
    }
}
