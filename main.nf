#!/usr/bin/env nextflow
nextflow.enable.dsl=2


//-- Help Message ---------------------------------------------------------------

def helpMessage() {
  log.info """
  ===============================
  nf-sra_screen
  Nextflow pipeline for screening SRA genomes
  Version: ${params.version}
  Author : Akito Shima (ASUQ)
  Email: asuq.4096@gmail.com
  ===============================
  Usage: nextflow run main.nf [parameters]

  Required parameters (assembly + binning):
    (A) SRA mode:
    --sra           Path to sra.csv (header: sra)
    (B) FASTQ mode:
    --fastq_tsv     Path to fastq.tsv (header: sample,read_type,reads)
                    (reads: comma-separated FASTQ paths)
    You can combine both modes by providing both --sra and --fastq_tsv

    And:
    --taxdump       Path to taxdump database folder
    --uniprot_db    Path to Uniprot database (.dmnd)

  Optional parameters (enable target taxa screening & extraction):
    --taxa          Path to taxa.csv for extraction (header: rank, taxa)
    --gtdb_ncbi_map Path to folder with GTDB-NCBI mapping Excel files
    --sandpiper_db  Path to Sandpiper database folder
    --singlem_db    Path to SingleM database folder

  Misc:
    --help          Show this help message
    --outdir        Output directory (default: ./output)
    --max_retries   Maximum number of retries for each process (default: 3)
  """.stripIndent()
}


def missingParametersError() {
    log.error "Missing input parameters"
    helpMessage()
    error """
    For assembly + binning, please provide:
      --taxdump and --uniprot_db
      and at least one of --sra or --fastq_tsv
    If you also want to enable target taxa screening & extraction, please provide:
      --taxa, --gtdb_ncbi_map, and --singlem_db
      and for sra mode also --sandpiper_db
    """.stripIndent()
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


// process COMEBIN {
//     tag "${sra}:${srr}"
//     label 'binning'
//     publishDir "${params.outdir}/${sra}/${srr}/binning",
//       mode: 'copy',
//       overwrite: true

//     input:
//     tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path(assembly_fasta), path(assembly_bam), path(assembly_csi)

//     output:
//     tuple val(sra), val(srr), path("comebin"),                                                                emit: bins
//     tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("comebin.note"), emit: note

//     script:
//     """
//     run_comebin_nf.sh \\
//       --assembly "${assembly_fasta}" \\
//       --bam "${assembly_bam}" \\
//       --cpus ${task.cpus} \\
//       --attempt ${task.attempt} \\
//       --max-retries ${params.max_retries}
//     """
// }


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
          path(semibin_dir),
          path(semibin_contig_bins),
          path(rosella_dir)

    output:
    tuple val(sra), val(srr), path("dastool"),                                                                emit: bins
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("dastool.note"), emit: note

    script:
    def metaDir     = metabat_dir          ?: ''
    // def comebinDir  = comebin_dir          ?: ''
    def semibinDir  = semibin_dir          ?: ''
    def semibinMap  = semibin_contig_bins  ?: ''
    def rosellaDir  = rosella_dir          ?: ''
    """
    run_dastool.sh \\
      --assembly "${assembly_fasta}" \\
      --metabat-dir "${metaDir}" \\
      --semibin-dir "${semibinDir}" \\
      --semibin-map "${semibinMap}" \\
      --rosella-dir "${rosellaDir}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries}
    """
}


process CREATE_EMPTY_SUMMARY {
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
    def doScreening = params.taxa != null

    // Start PRE_SCREENING after taxa validation
    def gated_sra_ch = doScreening
                        ? sra_ch.combine(validated_taxa).map { sra, vt -> sra }
                        : sra_ch

    // Extract metadata & filter SRR
    sra_metadata = DOWNLOAD_SRA_METADATA(gated_sra_ch)
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

    def singlem_reads_ch
    def sandpiper_note_ch
    def download_srr_note_ch
    def singlem_note_ch

    if (doScreening) {
      // SANDPIPER prescreening
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

      // Download SRR reads
      download_srr = DOWNLOAD_SRR(srr_prescreened)

      srr_reads = download_srr.reads.map { sra, srr, platform, model, strategy, asm, sandpiper_dec, reads, asm_txt ->
        def fixedAsm = file(asm_txt)?.text?.trim() ?: asm
        tuple(sra, srr, platform, model, strategy, fixedAsm, sandpiper_dec, reads)
      }

      // SINGLEM prescreening
      singlem = SINGLEM(srr_reads, validated_taxa, singlem_db_ch)
      singlem_reads_ch     = singlem.reads

      sandpiper_note_ch    = sandpiper.note
      download_srr_note_ch = download_srr.note
      singlem_note_ch      = singlem.note
    }

    else {
      // If no screening, set all decisions to PASS
      srr_prescreened = srr_ch.map { sra, srr, platform, model, strategy, assembler ->
        tuple(sra, srr, platform, model, strategy, assembler, 'PASS')
      }

      // Download SRR reads
      download_srr = DOWNLOAD_SRR(srr_prescreened)

      // normalise assembler, then drop sandpiper decision
      singlem_reads_ch = download_srr.reads.map { sra, srr, platform, model, strategy, asm, sandpiper_dec, reads, asm_txt ->
        def fixedAsm = file(asm_txt)?.text?.trim() ?: asm
        tuple(sra, srr, platform, model, strategy, fixedAsm, reads)
      }

      sandpiper_note_ch    = channel.empty()
      download_srr_note_ch = download_srr.note
      singlem_note_ch      = channel.empty()
    }

  emit:
    // For assembly
    singlem_reads         = singlem_reads_ch

    // For summary
    sra_metadata_skipped  = sra_metadata.skipped_sra
    sra_metadata_note     = sra_metadata.note
    sandpiper_note        = sandpiper_note_ch
    download_srr_note     = download_srr_note_ch
    singlem_note          = singlem_note_ch
}


workflow ASSEMBLY {
  take:
    singlem_reads
    validated_taxa
    uniprot_db_ch
    taxdump_ch

  main:
    def doScreening = params.taxa != null

    // Step 5: assemble reads
    short_ch    = singlem_reads.filter { sra, srr, platform, model, strategy, assembler, reads -> assembler.equalsIgnoreCase('short') }
    nanopore_ch = singlem_reads.filter { sra, srr, platform, model, strategy, assembler, reads -> assembler.equalsIgnoreCase('nanopore') }
    pacbio_ch   = singlem_reads.filter { sra, srr, platform, model, strategy, assembler, reads -> assembler.equalsIgnoreCase('pacbio') }
    hifi_ch     = singlem_reads.filter { sra, srr, platform, model, strategy, assembler, reads -> assembler.equalsIgnoreCase('hifi') }

    spades_asm     = METASPADES(short_ch)
    flyenano_asm   = METAFLYE_NANO(nanopore_ch)
    flyepacbio_asm = METAFLYE_PACBIO(pacbio_ch)
    myloasm_asm    = MYLOASM(hifi_ch)

    // Step 6: DIAMOND
    asm_fasta_ch = channel.empty()
                        .mix(spades_asm.assembly_fasta)
                        .mix(flyenano_asm.assembly_fasta)
                        .mix(flyepacbio_asm.assembly_fasta)
                        .mix(myloasm_asm.assembly_fasta)

    diamond = DIAMOND(asm_fasta_ch, uniprot_db_ch)

    // BAMs for BlobTools and binning
    bam_src = channel.empty()
                  .mix(spades_asm.assembly_bam)
                  .mix(flyenano_asm.assembly_bam)
                  .mix(flyepacbio_asm.assembly_bam)
                  .mix(myloasm_asm.assembly_bam)

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
    def taxa_summary_ch
    def taxa_note_ch

    if (doScreening) {
      taxa_extraction = EXTRACT_TAXA(blobtools.blobtable, validated_taxa, taxdump_ch)
      taxa_summary_ch = taxa_extraction.summary
      taxa_note_ch    = taxa_extraction.note

    }
    else{
      // No-taxa mode:
      // Treat any sample that reached BlobTools as "successful" and
      // synthesize an empty per-sample summary.csv using CREATE_EMPTY_SUMMARY.

      no_taxa_success_meta = blobtools.blobtable
        .map { sra, srr, platform, model, strategy, assembler, assembly_fasta, blobtable_csv ->
          // note field is empty string for successful runs
          tuple(sra, srr, platform, model, strategy, assembler, '')
        }
      no_taxa_summary = CREATE_EMPTY_SUMMARY(no_taxa_success_meta).skipped_rows

      taxa_summary_ch = no_taxa_summary.map { sra, srr, platform, model, strategy, assembler, summary_csv, note ->
        tuple(sra, srr, platform, model, strategy, assembler, summary_csv)
      }

      taxa_note_ch    = channel.empty()
    }

    assembly_notes_ch = channel.empty()
                          .mix(spades_asm.note)
                          .mix(flyenano_asm.note)
                          .mix(flyepacbio_asm.note)
                          .mix(myloasm_asm.note)

  emit:
    // For binning
    assembly_bam_all = bam_src
    blobtable        = blobtools.blobtable

    // For summary
    taxa_summary     = taxa_summary_ch
    taxa_note        = taxa_note_ch
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
    // comebin_binning = COMEBIN(binning_input)
    semibin_binning = SEMIBIN(binning_input, uniprot_db_ch)
    rosella_binning = ROSELLA(binning_input)

    // Prepare DASTool input
    dastool_base_by = binning_input.map { sra, srr, platform, model, strategy, assembler, assembly_fasta, bam, csi ->
      tuple([sra, srr], [platform, model, strategy, assembler, assembly_fasta])
    }

    // Tag each entry with a label
    metabat_entries = metabat_binning.bins.map { sra, srr, metabat_dir ->
      tuple([sra, srr], metabat_dir)
    }
    // comebin_entries = comebin_binning.bins.map { sra, srr, comebin_dir ->
    //   tuple([sra, srr], comebin_dir)
    // }
    semibin_entries = semibin_binning.bins.map { sra, srr, semibin_dir, semibin_bins ->
      tuple([sra, srr], [semibin_dir, semibin_bins])
    }
    rosella_entries = rosella_binning.bins.map { sra, srr, rosella_dir ->
      tuple([sra, srr], rosella_dir)
    }

    // Group all tool-entries by sample
    dastool_join = dastool_base_by
      .join(metabat_entries)
      .join(semibin_entries)
      .join(rosella_entries)
      // .join(comebin_entries)

    // Build DASTool input
    dastool_in = dastool_join.map { key, meta, metabat_dir, semibin_pair, rosella_dir ->
      def (sra, srr) = key
      def (platform, model, strategy, assembler, assembly_fasta) = meta
      def (semibin_dir, semibin_bins) = semibin_pair
      tuple(
        sra, srr, platform, model, strategy, assembler, assembly_fasta,
        metabat_dir, semibin_dir, semibin_bins, rosella_dir
      )
    }
    .filter { it != null }

    dastool_binning = DASTOOL(dastool_in)

  emit:
    // Expose raw note channels so SUMMARY can aggregate them together
    metabat_note = metabat_binning.note
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

    // from BINNING (empty channels when binning disabled)
    metabat_note
    semibin_note
    rosella_note
    dastool_note

    // constant
    outdir

  main:
    def doScreening = (params.taxa != null)
    def doBinning   = (params.binning == true)

    // 1) Convert skipped SRA metadata into per‑SRR rows with a textual note
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


    // 2) Base "successful" samples: anything with a taxa_summary (summary.csv)
    //    Note starts empty here; will add taxa/binning notes later.
    succeeded_sra = taxa_summary.map { sra, srr, platform, model, strategy, assembler, summary_csv ->
      tuple(sra, srr, platform, model, strategy, assembler, summary_csv, '')
    }


    // 3) Attach taxa notes and classify soft vs fatal EXTRACT_TAXA outcomes
    def succeeded_with_taxa = succeeded_sra
    def taxa_fatal_errors   = channel.empty()

    if( doScreening ) {
      // Exactly 2 entries per sample:
      // - [summary_csv, base_note] from succeeded_sra
      // - note_path from taxa_note
      def N_TAXA_ITEMS = 2

      succ_keyed = succeeded_sra.map { sra, srr, platform, model, strategy, assembler, summary_csv, base_note ->
        def keyMap = [sra: sra, srr: srr, platform: platform,
                      model: model, strategy: strategy, assembler: assembler]
        def key = groupKey(keyMap, N_TAXA_ITEMS)
        tuple(key, [summary_csv, base_note])
      }

      taxa_keyed = taxa_note.map { sra, srr, platform, model, strategy, assembler, note_path ->
        def keyMap = [sra: sra, srr: srr, platform: platform,
                      model: model, strategy: strategy, assembler: assembler]
        def key = groupKey(keyMap, N_TAXA_ITEMS)
        tuple(key, note_path)
      }

      succ_and_taxa = channel.empty()
        .mix(succ_keyed)
        .mix(taxa_keyed)
        .groupTuple()

      // Unpack to: (sra, srr, platform, model, strategy, assembler, summary_csv, base_note, taxa_text)
      succ_and_taxa_annot = succ_and_taxa
        .map { key, values ->
          def m = (Map) key.target
          def (sra, srr, platform, model, strategy, assembler) =
            [m.sra, m.srr, m.platform, m.model, m.strategy, m.assembler]

          // [summary_csv, base_note]
          def summaryEntry = values.find { it instanceof List && it.size() == 2 }
          if( !summaryEntry ) {
            return null
          }

          def summary_csv = summaryEntry[0]
          def base_note   = summaryEntry[1] ?: ''

          // note_path for EXTRACT_TAXA
          def taxa_note_file = values.find { !(it instanceof List) }
          def taxa_text = ''
          if( taxa_note_file ) {
            taxa_text = file(taxa_note_file as String).text.trim()
          }

          tuple(sra, srr, platform, model, strategy, assembler, summary_csv, base_note, taxa_text)
        }
        .filter { it != null }

      // Split EXTRACT_TAXA outcomes into "success" vs "fatal"
      def softPattern = 'skipping extraction because taxa list contains only GTDB-style taxa'

      def taxa_branches = succ_and_taxa_annot.branch { sra, srr, platform, model, strategy, assembler, summary_csv, base_note, taxa_text ->
        // Empty note or GTDB-only “soft skip” => success
        success: (!taxa_text || taxa_text.contains(softPattern))
        // Anything else => fatal
        fatal:   (taxa_text && !taxa_text.contains(softPattern))
      }

      taxa_success = taxa_branches.success
      taxa_fatal   = taxa_branches.fatal

      // Non-fatal EXTRACT_TAXA results remain as "successful" samples.
      // Append taxa_text (if any) to the note field.
      succeeded_with_taxa = taxa_success.map { sra, srr, platform, model, strategy, assembler, summary_csv, base_note, taxa_text ->
        def final_note = base_note
        if( taxa_text ) {
          final_note = final_note ? "${final_note}; ${taxa_text}" : taxa_text
        }
        tuple(sra, srr, platform, model, strategy, assembler, summary_csv, final_note)
      }

      // Fatal EXTRACT_TAXA outcomes become "errors" that will go through
      // CREATE_EMPTY_SUMMARY, similar to other fatal notes.
      taxa_fatal_errors = taxa_fatal.map { sra, srr, platform, model, strategy, assembler, summary_csv, base_note, taxa_text ->
        def note_text = base_note
        if( taxa_text ) {
          note_text = note_text ? "${note_text}; ${taxa_text}" : taxa_text
        }
        tuple(sra, srr, platform, model, strategy, assembler, note_text)
      }
    }


    // 4) Aggregate binning notes into a single string per sample (if binning)
    def binning_agg = channel.empty()

    if( doBinning ) {
      // Expect exactly 4 notes per sample: metabat, semibin, rosella, dastool.
      def N_BIN_STATUS = 4

      metabat_status = metabat_note.map { sra, srr, platform, model, strategy, assembler, note_path ->
        def keyMap = [sra: sra, srr: srr, platform: platform, model: model,
                      strategy: strategy, assembler: assembler]
        def key = groupKey(keyMap, N_BIN_STATUS)
        tuple(key, ['metabat', note_path])
      }

      semibin_status = semibin_note.map { sra, srr, platform, model, strategy, assembler, note_path ->
        def keyMap = [sra: sra, srr: srr, platform: platform, model: model,
                      strategy: strategy, assembler: assembler]
        def key = groupKey(keyMap, N_BIN_STATUS)
        tuple(key, ['semibin', note_path])
      }

      rosella_status = rosella_note.map { sra, srr, platform, model, strategy, assembler, note_path ->
        def keyMap = [sra: sra, srr: srr, platform: platform, model: model,
                      strategy: strategy, assembler: assembler]
        def key = groupKey(keyMap, N_BIN_STATUS)
        tuple(key, ['rosella', note_path])
      }

      dastool_status = dastool_note.map { sra, srr, platform, model, strategy, assembler, note_path ->
        def keyMap = [sra: sra, srr: srr, platform: platform, model: model,
                      strategy: strategy, assembler: assembler]
        def key = groupKey(keyMap, N_BIN_STATUS)
        tuple(key, ['dastool', note_path])
      }

      binning_mix = channel.empty()
        .mix(metabat_status)
        .mix(semibin_status)
        .mix(rosella_status)
        .mix(dastool_status)
        .groupTuple()

      // binning_agg: (sra, srr, platform, model, strategy, assembler, "tool1: msg; tool2: msg; ...")
      binning_agg = binning_mix.map { key, items ->
        def m = (Map) key.target
        def (sra, srr, platform, model, strategy, assembler) =
          [m.sra, m.srr, m.platform, m.model, m.strategy, m.assembler]

        def note_texts = items.collect { entry ->
          def (tool, note_path) = entry
          def txt = file(note_path as String).text.trim()
          txt ?: null
        }.findAll { it }

        def joined = note_texts ? note_texts.join('; ') : ''
        tuple(sra, srr, platform, model, strategy, assembler, joined)
      }
    }


    // 5) Combine successful rows with aggregated binning notes (if binning)
    def final_success = succeeded_with_taxa

    if( doBinning ) {
      def N_SUCCESS_ITEMS = 2

      succ_keyed2 = succeeded_with_taxa.map { sra, srr, platform, model, strategy, assembler, summary_csv, base_note ->
        def keyMap = [sra: sra, srr: srr, platform: platform,
                      model: model, strategy: strategy, assembler: assembler]
        def key = groupKey(keyMap, N_SUCCESS_ITEMS)
        tuple(key, [summary_csv, base_note])
      }

      binning_keyed = binning_agg.map { sra, srr, platform, model, strategy, assembler, binning_note ->
        def keyMap = [sra: sra, srr: srr, platform: platform,
                      model: model, strategy: strategy, assembler: assembler]
        def key = groupKey(keyMap, N_SUCCESS_ITEMS)
        tuple(key, binning_note)
      }

      succ_and_bin = channel.empty()
        .mix(succ_keyed2)
        .mix(binning_keyed)
        .groupTuple()

      final_success = succ_and_bin
        .map { key, values ->
          def m = (Map) key.target
          def (sra, srr, platform, model, strategy, assembler) =
            [m.sra, m.srr, m.platform, m.model, m.strategy, m.assembler]

          def summaryEntry = values.find { it instanceof List && it.size() == 2 }
          if( !summaryEntry ) {
            return null
          }

          def summary_csv = summaryEntry[0]
          def base_note   = summaryEntry[1] ?: ''

          def binning_note_text = ''
          def extra = values.find { !(it instanceof List) } as String
          if( extra ) {
            binning_note_text = extra.trim()
          }

          def final_note = base_note
          if( binning_note_text ) {
            final_note = final_note ? "${final_note}; ${binning_note_text}" : binning_note_text
          }

          tuple(sra, srr, platform, model, strategy, assembler, summary_csv, final_note)
        }
        .filter { it != null }
    }


    // 6) Collect all fatal errors (including fatal EXTRACT_TAXA) and turn them
    //    into empty per-sample summaries.
    errors = channel.empty()
      .mix(sra_metadata_note)
      .mix(sandpiper_note)
      .mix(download_srr_note)
      .mix(singlem_note)
      .mix(assembly_notes)
      .mix(diamond_note)
      .mix(blobtools_note)
      .map { sra, srr, platform, model, strategy, assembler, note_path ->
        def note = file(note_path).text.trim()
        tuple(sra, srr, platform, model, strategy, assembler, note)
      }
      // Fatal EXTRACT_TAXA errors (inc. "run failed") are added here
      .mix(taxa_fatal_errors)
      // Plus skipped SRR rows
      .mix(skipped_srr)

    // failed_sra: (sra, srr, platform, model, strategy, assembler, empty_summary.csv, note)
    failed_sra = CREATE_EMPTY_SUMMARY(errors)


    // 7) Combine succeeded and failed, and append rows to summary.tsv
    summary = channel.empty()
      .mix(final_success)
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

    def sraMode   = params.sra != null
    def fastqMode = params.fastq_tsv != null

    if (!sraMode && !fastqMode) {
      log.error "Error: either --sra or --fastq_tsv parameter must be provided"
      missingParametersError()
    }

    if (!params.taxdump || !params.uniprot_db) {
      log.error "Error: Missing --taxdump or --uniprot_db"
      missingParametersError()
    }

    def doScreening = params.taxa != null
    if (doScreening) {
      if (!params.taxdump || !params.gtdb_ncbi_map || !params.singlem_db) {
        log.error "Error: Missing --taxdump, --gtdb_ncbi_map, or --singlem_db required for taxa filtering"
        missingParametersError()
      }
      if (sraMode && !params.sandpiper_db) {
        log.error "Error: --sandpiper_db is required for SRA-based taxa screening"
        missingParametersError()
      }
    }

    def doBinning = (params.binning == true)

    // Common channels
    def outdir        = file(params.outdir).toAbsolutePath().toString()
    def taxdump_ch    = channel.value( file(params.taxdump) )
    def uniprot_db_ch = channel.value( file(params.uniprot_db) )

    // Taxa-related channels
    def validated_taxa_ch = channel.empty()
    def singlem_db_ch     = channel.empty()
    def sandpiper_db_ch   = channel.empty()
    if (doScreening) {
      def taxa_ch          = channel.value( file(params.taxa) )
      def gtdb_ncbi_map_ch = channel.value( file(params.gtdb_ncbi_map) )
      singlem_db_ch        = channel.value( file(params.singlem_db) )
      sandpiper_db_ch      = sraMode ? channel.value( file(params.sandpiper_db) ) : channel.empty()

      // validate taxa
      validated_taxa_ch = VALIDATE_TAXA(taxa_ch, taxdump_ch, gtdb_ncbi_map_ch).valid_taxa
    }


    // SRA mode
    def sra_singlem_reads_ch    = channel.empty()
    def sra_metadata_skipped_ch = channel.empty()
    def sra_metadata_note       = channel.empty()
    def sra_sandpiper_note      = channel.empty()
    def sra_download_srr_note   = channel.empty()
    def sra_singlem_note        = channel.empty()
    if (sraMode) {
      sra_ch = channel.fromPath(params.sra, checkIfExists: true)
                      .splitCsv(header: true, strip: true)
                      .map { row -> row.sra.trim() }
                      .filter { it }
                      .distinct()

      pre = PRE_SCREENING(sra_ch, validated_taxa_ch, sandpiper_db_ch, singlem_db_ch)

      sra_singlem_reads_ch    = pre.singlem_reads
      sra_metadata_skipped_ch = pre.sra_metadata_skipped
      sra_metadata_note       = pre.sra_metadata_note
      sra_sandpiper_note      = pre.sandpiper_note
      sra_download_srr_note   = pre.download_srr_note
      sra_singlem_note        = pre.singlem_note
    }

    // FASTQ mode
    def fastq_singlem_reads_ch = channel.empty()
    def fastq_singlem_note     = channel.empty()
    if (fastqMode) {
      fastq_ch = channel.fromPath(params.fastq_tsv, checkIfExists: true)
                        .splitCsv(header: true, sep: '\t', strip: true)
                        .map { row ->
                          def sample    = (row.sample ?: '').trim()
                          def read_type = (row.read_type ?: '').trim()
                          def reads_raw = (row.reads ?: '').trim()

                          if (!sample || !read_type || !reads_raw) {
                            log.warn "Skipping FASTQ TSV row with missing fields: ${row}"
                            return null
                          }

                          def read_files = reads_raw.split(/\s*,\s*/).findAll { it }.collect { file(it) }

                          if (!read_files) {
                            log.warn "No valid FASTQ paths for sample ${sample}; skipping"
                            return null
                          }

                          def sra       = sample
                          def srr       = sample
                          def platform  = "UNKNOWN"
                          def model     = read_type
                          def strategy  = "UNKNOWN"
                          def assembler = read_type

                          tuple(sra, srr, platform, model, strategy, assembler, read_files)
                        }
                        .filter { it != null }

      if (doScreening) {
        fastq_for_singlem_ch = fastq_ch.map { sra, srr, platform, model, strategy, assembler, reads ->
          tuple(sra, srr, platform, model, strategy, assembler, "RUN_SINGLEM", reads)
        }

        singlem_fastq = SINGLEM(fastq_for_singlem_ch, validated_taxa_ch, singlem_db_ch)
        fastq_singlem_reads_ch = singlem_fastq.reads
        fastq_singlem_note     = singlem_fastq.note
      }
      else {
        fastq_singlem_reads_ch = fastq_ch
        fastq_singlem_note     = channel.empty()
      }
    }

    // Merge SRA and FASTQ reads for assembly
    def singlem_reads_all = channel.empty()
                                  .mix(sra_singlem_reads_ch)
                                  .mix(fastq_singlem_reads_ch)
    def singlem_notes_all = channel.empty()
                                  .mix(sra_singlem_note)
                                  .mix(fastq_singlem_note)

    asm = ASSEMBLY(singlem_reads_all, validated_taxa_ch, uniprot_db_ch, taxdump_ch)

    // BINNING notes: default to empty channels when binning disabled
    def metabat_note_ch = channel.empty()
    def semibin_note_ch = channel.empty()
    def rosella_note_ch = channel.empty()
    def dastool_note_ch = channel.empty()

    if (doBinning) {
      def binning = BINNING(asm.blobtable, asm.assembly_bam_all, uniprot_db_ch)
      metabat_note_ch = binning.metabat_note
      semibin_note_ch = binning.semibin_note
      rosella_note_ch = binning.rosella_note
      dastool_note_ch = binning.dastool_note
    }

    SUMMARY(
      // PRE_SCREENING: summary-related outputs
      sra_metadata_skipped_ch,
      sra_metadata_note,
      sra_sandpiper_note,
      sra_download_srr_note,
      singlem_notes_all,

      // ASSEMBLY: summary-related outputs
      asm.assembly_notes,
      asm.diamond_note,
      asm.blobtools_note,
      asm.taxa_note,
      asm.taxa_summary,

      // BINNING: raw notes from each binner (or empty channels if binning disabled)
      metabat_note_ch,
      semibin_note_ch,
      rosella_note_ch,
      dastool_note_ch,

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
