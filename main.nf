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
    publishDir "${params.outdir}/${sra}/${srr}/binning",
      mode: 'copy',
      overwrite: true,
      saveAs: { filename ->
        if (filename == "FAIL.note") return "metabat.FAIL.note"
        return filename
      }

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path(assembly_fasta), path(assembly_bam)

    output:
    tuple val(sra), val(srr), path("metabat"),                                                             optional: true, emit: bins
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("FAIL.note"), optional: true, emit: note

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
    publishDir "${params.outdir}/${sra}/${srr}/binning",
      mode: 'copy',
      overwrite: true,
      saveAs: { filename ->
        if (filename == "FAIL.note") return "concoct.FAIL.note"
        return filename
      }

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path(assembly_fasta), path(assembly_bam)

    output:
    tuple val(sra), val(srr), path("concoct"),                                                             optional: true, emit: bins
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("FAIL.note"), optional: true, emit: note

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
    publishDir "${params.outdir}/${sra}/${srr}/binning",
      mode: 'copy',
      overwrite: true,
      saveAs: { filename ->
        if (filename == "FAIL.note") return "semibin.FAIL.note"
        return filename
      }

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path(assembly_fasta), path(assembly_bam)
    path uniprot_db

    output:
    tuple val(sra), val(srr), path("semibin"), path("semibin/contig_bins.tsv"),                            optional: true, emit: bins
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("FAIL.note"), optional: true, emit: note

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
    publishDir "${params.outdir}/${sra}/${srr}/binning",
      mode: 'copy',
      overwrite: true,
      saveAs: { filename ->
        if (filename == "FAIL.note") return "rosella.FAIL.note"
        return filename
      }

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path(assembly_fasta), path(assembly_bam)

    output:
    tuple val(sra), val(srr), path("rosella"),                                                             optional: true, emit: bins
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("FAIL.note"), optional: true, emit: note

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
      overwrite: true,
      saveAs: { filename ->
        if (filename == "FAIL.note") return "dastool.FAIL.note"
        return filename
      }

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler),
          path(assembly_fasta),
          val(metabat_dir),
          val(concoct_dir),
          val(semibin_dir),
          val(semibin_contig_bins),
          val(rosella_dir)

    output:
    tuple val(sra), val(srr), path("dastool"),                                                             optional: true, emit: bins
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("FAIL.note"), optional: true, emit: note

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
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path(note_files)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("binning_note.txt"), emit: binning_notes

    script:
    // note_files is a List<File> staged into the work dir
    def list = note_files.collect { it.name }.join(' ')
    """
    # Concatenate and flatten all binning FAIL.note messages into a single line
    cat ${list} 2>/dev/null \\
      | tr '\\n' ' ' \\
      | sed 's/  */ /g' \\
      > binning_note.txt
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
    bam_by = assembly_bam_ch.map { sra, srr, bam, csi -> tuple([sra, srr], bam) }
    binning_base_by = blobtable_ch.map { sra, srr, platform, model, strategy, assembler, assembly_fasta, blobtable ->
      tuple([sra, srr], [platform, model, strategy, assembler, assembly_fasta])
    }
    binning_join = binning_base_by.join(bam_by)

    // binning_input: (sra, srr, platform, model, strategy, assembler, assembly_fasta, bam)
    binning_input = binning_join.map { key, meta, bam ->
      def (sra, srr) = key
      def (platform, model, strategy, assembler, assembly_fasta) = meta
      tuple(sra, srr, platform, model, strategy, assembler, assembly_fasta, bam)
    }

    // Run individual binners
    metabat_binning = METABAT(binning_input)
    concoct_binning = CONCOCT(binning_input)
    semibin_binning = SEMIBIN(binning_input, uniprot_db_ch)
    rosella_binning = ROSELLA(binning_input)

    // Prepare DASTool input: run if at least ONE binner produced bins
    dastool_base_by = binning_input.map { sra, srr, platform, model, strategy, assembler, assembly_fasta, bam ->
      tuple([sra,srr], [platform, model, strategy, assembler, assembly_fasta])
    }

    // Tag each entry with a label
    meta_by      = dastool_base_by.map { key, meta -> tuple(key, ['meta', meta]) }
    metabat_by_key = metabat_binning.bins.map { sra, srr, metabat_dir -> tuple([sra,srr], ['metabat', metabat_dir]) }
    concoct_by_key = concoct_binning.bins.map { sra, srr, concoct_dir -> tuple([sra,srr], ['concoct', concoct_dir]) }
    semibin_by_key = semibin_binning.bins.map { sra, srr, semibin_dir, contig_bins -> tuple([sra,srr], ['semibin', [semibin_dir, contig_bins]]) }
    rosella_by_key = rosella_binning.bins.map { sra, srr, rosella_dir -> tuple([sra,srr], ['rosella', rosella_dir]) }

    // Group all tool-entries by sample
    grouped = channel.empty()
                .mix(meta_by)
                .mix(metabat_by_key)
                .mix(concoct_by_key)
                .mix(semibin_by_key)
                .mix(rosella_by_key)
                .groupTuple()

    // Build DASTool input:
    // - keep only samples that have at least one tool (metabat/concoct/semibin/rosella)
    // - allow others to be null
    dastool_in = grouped.map { key, entries ->
        def (sra, srr) = key
        def meta_entry = entries.find { it[0] == 'meta' }
        if (!meta_entry) return null  // should not happen

        def meta = meta_entry[1]
        def (platform, model, strategy, assembler, assembly_fasta) = meta
        def metabat_dir  = entries.find { it[0] == 'metabat' }?.getAt(1)
        def concoct_dir  = entries.find { it[0] == 'concoct' }?.getAt(1)
        def semibin_data = entries.find { it[0] == 'semibin' }?.getAt(1)
        def rosella_dir  = entries.find { it[0] == 'rosella' }?.getAt(1)

        // If none of the tools produced bins, skip this sample
        if (!metabat_dir && !concoct_dir && !semibin_data && !rosella_dir)
          return null

        def (semibin_dir, contig_bins) = semibin_data ?: [null, null]

        tuple(
          sra, srr, platform, model, strategy, assembler,
          assembly_fasta,
          metabat_dir,
          concoct_dir,
          semibin_dir,
          contig_bins,
          rosella_dir
        )
      }
      .filter { it != null }

    dastool_binning = DASTOOL(dastool_in)

  emit:
    metabat_note = metabat_binning.note
    concoct_note = concoct_binning.note
    semibin_note = semibin_binning.note
    rosella_note = rosella_binning.note
    dastool_note = dastool_binning.note
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
    metabat_note
    concoct_note
    semibin_note
    rosella_note
    dastool_note

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

    // Collect all binning FAIL.notes (MetaBAT, CONCOCT, SemiBin2, Rosella, DASTool)
    binning_errors_raw = channel.empty()
                          .mix(metabat_note)
                          .mix(concoct_note)
                          .mix(semibin_note)
                          .mix(rosella_note)
                          .mix(dastool_note)

    // Group binning FAIL.notes by sample
    binning_errors_grouped = binning_errors_raw
      .map { sra, srr, platform, model, strategy, assembler, note_path ->
        def key = [sra, srr, platform, model, strategy, assembler]
        tuple(key, note_path)
      }
      .groupTuple()
      .map { key, note_paths ->
        def (sra, srr, platform, model, strategy, assembler) = key
        tuple(sra, srr, platform, model, strategy, assembler, note_paths)
      }

    // Aggregate binning errors per sample into a single note string
    binning_notes = BINNING_ERROR_SUMMARY(binning_errors_grouped).binning_notes

    // Combine succeeded_sra and binning_notes by key
    succ_keyed = succeeded_sra.map { sra, srr, platform, model, strategy, assembler, summary_csv, note ->
      def key   = [sra, srr, platform, model, strategy, assembler]
      def value = [summary_csv, note]  // base note is ''
      tuple(key, value)
    }
    binning_keyed = binning_notes.map { sra, srr, platform, model, strategy, assembler, bin_note_file ->
      def key  = [sra, srr, platform, model, strategy, assembler]
      def note = file(bin_note_file).text.trim()
      tuple(key, note)
    }
    succ_and_bin = channel.empty()
                        .mix(succ_keyed)
                        .mix(binning_keyed)
                        .groupTuple()

    // Build final succeeded_sra with binning annotations in the note field
    succeeded_with_binning = succ_and_bin
      .map { key, values ->
        def (sra, srr, platform, model, strategy, assembler) = key

        // Find the summary entry (List [summary_csv, base_note])
        def summaryEntry = values.find { it instanceof List && it.size() == 2 }
        if (!summaryEntry) {
          return null   // no summary => not a "successful" sample; ignore here
        }

        def summary_csv = summaryEntry[0]
        def base_note   = summaryEntry[1] ?: ''

        // All other values are binning notes (Strings)
        def bin_notes = values.findAll { !(it instanceof List) }
                              .collect { it as String }
                              .findAll { it }

        def bin_note   = bin_notes ? bin_notes.join('; ') : ''
        def final_note = base_note
        if (bin_note) {
          final_note = final_note ? "${final_note}; ${bin_note}" : bin_note
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
  main:
    if (params.help) {
      helpMessage()
      exit 0
    }

    if (!params.sra || !params.uniprot_db || !params.taxa \
        || !params.taxdump || !params.gtdb_ncbi_map \
        || !params.sandpiper_db || !params.singlem_db) {
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
      binning.metabat_note,
      binning.concoct_note,
      binning.semibin_note,
      binning.rosella_note,
      binning.dastool_note,

      // constant
      outdir
    )


  onComplete:
    def outdirPath  = file(params.outdir ?: './output').toAbsolutePath()
    def summaryFile = outdirPath.resolve('summary.tsv')
    def traceFile   = file("${workflow.launchDir}/execution-report/trace.tsv").toAbsolutePath()
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
