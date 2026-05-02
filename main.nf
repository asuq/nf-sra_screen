#!/usr/bin/env nextflow

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
    --noassembly    Skip ASSEMBLY/BINNING
    --binning       Also run BINNING after ASSEMBLY
    --assemblers    Assembly tools: auto, all, or comma-separated tool names
                    (default: auto; aliases: --assembler, spades, flye)
    --binners       Comma-separated binners (default: all-compatible; allowed: all-compatible,auto,metabat,semibin,rosella,comebin,vamb,lorbin)
    --refiners      Comma-separated refiners (default: binette; allowed: dastool,binette)
    --checkm2_db    CheckM2 DIAMOND database required with --refiners binette
    --semibin_environment  SemiBin2 pretrained environment (default: global)
    --gpu           Use GPU variants for COMEBin, VAMB, and HiFi-only LorBin
    --gpu_type      Optional GPU type for typed scheduler requests on GWDG
    --gpus          GPU count for scheduler requests on GWDG (default: 1)
    --outdir        Output directory (default: ./output)
    --max_retries   Maximum number of retries for each process (default: 3)
  """.stripIndent()
}


def missingParametersError() {
  log.error "Missing input parameters"
  helpMessage()

  def noAssembly = params.noassembly?.toString()?.toBoolean() ?: false

  if (noAssembly) {
    error """
    For --noassembly (only PRE_SCREENING), please provide:
      at least one of --sra or --fastq_tsv
    If you enable taxa screening (--taxa), please also provide:
      --taxdump, --gtdb_ncbi_map, and --singlem_db
      and for sra mode also --sandpiper_db
    """.stripIndent()
  }
  else {
    error """
    For assembly + binning, please provide:
      --taxdump and --uniprot_db
      and at least one of --sra or --fastq_tsv
    If you also want to enable target taxa screening & extraction, please provide:
      --taxa, --gtdb_ncbi_map, and --singlem_db
      and for sra mode also --sandpiper_db
    """.stripIndent()
  }
}


include {
  validateBinningOptions
  binnersForSample
  binnerCsvContains
  binnerCsvSize
  selectedAssemblerTokens
  assemblersForReadType
  assemblerPublishDir
} from './lib/workflow_helpers.nf'


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
    publishDir { "${params.outdir}/metadata/${sra}/" },
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
    tuple val(sra), path("${sra}.filtered.csv"),                                             optional: true, emit: filtered_sra
    tuple val(sra), path("${sra}.skipped.csv"),                                              optional: true, emit: skipped_sra
    tuple val(sra), val(''), val(''), val(''), val(''), val(''), val(''), path("FAIL.note"), optional: true, emit: note

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
    publishDir { "${params.outdir}/${sra}/${srr}/" }, mode: 'copy', overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type)
    path valid_taxa
    path sandpiper_db

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), path("sandpiper_decision.txt"),            emit: decision
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(''), path("FAIL.note"), optional: true, emit: note
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
    publishDir { "${params.outdir}/${sra}/${srr}/" },
      mode: 'copy',
      overwrite: true,
      saveAs: { filename ->
        filename in ["FAIL.note"] ? filename : null
      }

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(sandpiper)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(sandpiper), path("*.f*q*"), path("assembler.txt"), optional: true, emit: reads
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(''), path("FAIL.note"),                            optional: true, emit: note

    script:
    """
    run_download_srr.sh \\
      --srr "${srr}" \\
      --platform "${platform}" \\
      --read-type "${read_type}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries}
    """
}


process SINGLEM {
    tag "${sra}:${srr}"
    publishDir { "${params.outdir}/${sra}/${srr}/" },
      mode: 'copy',
      overwrite: true,
      saveAs: { filename ->
        if (filename == "singlem_output.tsv") return filename
        if (filename.startsWith("singlem_taxonomic_profile")) return filename
        if (filename == "FAIL.note") return filename
        return null
      }

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(sandpiper), path(reads)
    path valid_taxa
    path singlem_db

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), path("reads_ok/*.f*q*"),          optional: true, emit: reads
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(''), path("FAIL.note"),       optional: true, emit: note
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


process FASTP_SHORT {
    tag "${sra}:${srr}:${assembler}"
    publishDir { assemblerPublishDir(sra, srr, assembler) },
      mode: 'copy',
      overwrite: true,
      saveAs: { filename ->
        filename in [
          "fastp.html",
          "fastp.json",
          "FAIL.note"
        ] ? filename : null
      }

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(reads)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("*_fastp_R*.fastq.gz"), optional: true, emit: reads
    tuple val(sra), val(srr), val(read_type), val(assembler), path("fastp.html"),                                                        optional: true, emit: html
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("FAIL.note"),               optional: true, emit: note

    script:
    """
    run_fastp_short.sh \\
      --srr "${srr}" \\
      --assembler "${assembler}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries} \\
      --reads ${reads}
    """
}


process METASPADES {
    tag "${sra}:${srr}:${assembler}"
    label 'assembly'
    publishDir { assemblerPublishDir(sra, srr, assembler) },
      mode: 'copy',
      overwrite: true,
      saveAs: { filename ->
        filename in [
          "assembly.fasta",
          "assembly.gfa",
          "spades.log",
          "FAIL.note"
        ] ? filename : null
      }

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(reads)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("assembly.fasta"), optional: true, emit: assembly_fasta
    tuple val(sra), val(srr), val(read_type), val(assembler), path("assembly.gfa"),                                            optional: true, emit: assembly_graph
    tuple val(sra), val(srr), val(read_type), val(assembler), path("spades.log"),                                              optional: true, emit: assembly_log
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("FAIL.note"),     optional: true, emit: note

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


process UNICYCLER {
    tag "${sra}:${srr}:${assembler}"
    label 'assembly'
    publishDir { assemblerPublishDir(sra, srr, assembler) },
      mode: 'copy',
      overwrite: true,
      saveAs: { filename ->
        filename in [
          "assembly.fasta",
          "assembly.gfa",
          "spades.log",
          "FAIL.note"
        ] ? filename : null
      }

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(reads)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("assembly.fasta"), optional: true, emit: assembly_fasta
    tuple val(sra), val(srr), val(read_type), val(assembler), path("assembly.gfa"),                                            optional: true, emit: assembly_graph
    tuple val(sra), val(srr), val(read_type), val(assembler), path("spades.log"),                                              optional: true, emit: assembly_log
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("FAIL.note"),     optional: true, emit: note

    script:
    """
    # --reads should be the last argument and unquoted to capture all read files
    run_unicycler.sh \\
      --srr "${srr}" \\
      --cpus ${task.cpus} \\
      --memory-gb ${task.memory.toGiga()} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries} \\
      --reads ${reads}
    """
}


process METAFLYE_NANO {
    tag "${sra}:${srr}:${assembler}"
    label 'assembly'
    publishDir { assemblerPublishDir(sra, srr, assembler) },
      mode: 'copy',
      overwrite: true,
      saveAs: { filename ->
        filename in [
          "assembly.fasta",
          "assembly.gfa",
          "flye.log",
          "FAIL.note"
        ] ? filename : null
      }

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(reads)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("assembly.fasta"), optional: true, emit: assembly_fasta
    tuple val(sra), val(srr), val(read_type), val(assembler), path("assembly.gfa"),                                            optional: true, emit: assembly_graph
    tuple val(sra), val(srr), val(read_type), val(assembler), path("flye.log"),                                                optional: true, emit: assembly_log
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("FAIL.note"),     optional: true, emit: note

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
    tag "${sra}:${srr}:${assembler}"
    label 'assembly'
    publishDir { assemblerPublishDir(sra, srr, assembler) },
      mode: 'copy',
      overwrite: true,
      saveAs: { filename ->
        filename in [
          "assembly.fasta",
          "assembly.gfa",
          "flye.log",
          "FAIL.note"
        ] ? filename : null
      }

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(reads)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("assembly.fasta"), optional: true, emit: assembly_fasta
    tuple val(sra), val(srr), val(read_type), val(assembler), path("assembly.gfa"),                                            optional: true, emit: assembly_graph
    tuple val(sra), val(srr), val(read_type), val(assembler), path("flye.log"),                                                optional: true, emit: assembly_log
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("FAIL.note"),     optional: true, emit: note

    script:
    """
    run_metaflye_pacbio.sh \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries} \\
      --reads ${reads}
    """
}


process METAFLYE_HIFI {
    tag "${sra}:${srr}:${assembler}"
    label 'assembly'
    publishDir { assemblerPublishDir(sra, srr, assembler) },
      mode: 'copy',
      overwrite: true,
      saveAs: { filename ->
        filename in [
          "assembly.fasta",
          "assembly.gfa",
          "flye.log",
          "FAIL.note"
        ] ? filename : null
      }

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(reads)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("assembly.fasta"), optional: true, emit: assembly_fasta
    tuple val(sra), val(srr), val(read_type), val(assembler), path("assembly.gfa"),                                            optional: true, emit: assembly_graph
    tuple val(sra), val(srr), val(read_type), val(assembler), path("flye.log"),                                                optional: true, emit: assembly_log
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("FAIL.note"),     optional: true, emit: note

    script:
    """
    run_metaflye_hifi.sh \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries} \\
      --reads ${reads}
    """
}


process MYLOASM {
    tag "${sra}:${srr}:${assembler}"
    label 'assembly'
    publishDir { assemblerPublishDir(sra, srr, assembler) },
      mode: 'copy',
      overwrite: true,
      saveAs: { filename ->
        if (filename == "assembly.fasta") return filename
        if (filename == "assembly.gfa") return filename
        if (filename.startsWith("myloasm_") && filename.endsWith(".log")) return filename
        if (filename == "FAIL.note") return filename
        return null
      }

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(reads)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("assembly.fasta"), optional: true, emit: assembly_fasta
    tuple val(sra), val(srr), val(read_type), val(assembler), path("assembly.gfa"),                                            optional: true, emit: assembly_graph
    tuple val(sra), val(srr), val(read_type), val(assembler), path("myloasm_*.log"),                                           optional: true, emit: assembly_log
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("FAIL.note"),     optional: true, emit: note

    script:
    """
    run_myloasm_hifi.sh \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries} \\
      --reads ${reads}
    """
}


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
}


process DIAMOND {
    tag "${sra}:${srr}:${assembler}"
    publishDir { assemblerPublishDir(sra, srr, assembler) }, mode: 'copy', overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(assembly_fasta)
    path uniprot_db

    output:
    tuple val(sra), val(srr), val(read_type), val(assembler), path("assembly_vs_uniprot.tsv"),                                     optional: true, emit: blast
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("FAIL.note"), optional: true, emit: note

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
    tag "${sra}:${srr}:${assembler}"
    publishDir { assemblerPublishDir(sra, srr, assembler) }, mode: 'copy', overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(assembly_fasta), path(blast), path(assembly_bam), path(assembly_csi)
    path taxdump

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(assembly_fasta), path("blobtools.csv"), optional: true, emit: blobtable
    tuple val(sra), val(srr), val(read_type), val(assembler), path("blobtools*.svg"),                                                               optional: true, emit: blobplots
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("FAIL.note"),                         optional: true, emit: note

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
    tag "${sra}:${srr}:${assembler}"
    publishDir { assemblerPublishDir(sra, srr, assembler) }, mode: 'copy', overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(assembly_fasta), path(blobtable)
    path valid_taxa
    path taxdump

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("summary.csv"), optional: true, emit: summary
    tuple val(sra), val(srr), val(read_type), val(assembler), path("*.ids.csv"),                                            optional: true, emit: extracted_ids
    tuple val(sra), val(srr), val(read_type), val(assembler), path("*.fasta"),                                              optional: true, emit: extracted_fasta
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("FAIL.note"),  optional: true, emit: note

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
    tag "${sra}:${srr}:${assembler}"
    label 'binning'
    publishDir { "${assemblerPublishDir(sra, srr, assembler)}/binning" },
      mode: 'copy',
      overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(assembly_fasta), path(assembly_bam), path(assembly_csi)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler),
          val("metabat"), path("metabat"), path("metabat.contig2bin.tsv"), path("metabat.note"),             emit: result

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


process COMEBIN {
    tag "${sra}:${srr}:${assembler}"
    label 'binning'
    publishDir { "${assemblerPublishDir(sra, srr, assembler)}/binning" },
      mode: 'copy',
      overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(assembly_fasta), path(assembly_bam), path(assembly_csi)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler),
          val("comebin"), path("comebin"), path("comebin.contig2bin.tsv"), path("comebin.note"),             emit: result

    script:
    """
    run_comebin_nf.sh \\
      --assembly "${assembly_fasta}" \\
      --bam "${assembly_bam}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries}
    """
}


process COMEBIN_GPU {
    tag "${sra}:${srr}:${assembler}"
    label 'gpu'
    publishDir { "${assemblerPublishDir(sra, srr, assembler)}/binning" },
      mode: 'copy',
      overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(assembly_fasta), path(assembly_bam), path(assembly_csi)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler),
          val("comebin"), path("comebin"), path("comebin.contig2bin.tsv"), path("comebin.note"),             emit: result

    script:
    """
    run_comebin_nf.sh \\
      --assembly "${assembly_fasta}" \\
      --bam "${assembly_bam}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries} \\
      --require-cuda
    """
}


process VAMB {
    tag "${sra}:${srr}:${assembler}"
    label 'binning'
    publishDir { "${assemblerPublishDir(sra, srr, assembler)}/binning" },
      mode: 'copy',
      overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(assembly_fasta), path(assembly_bam), path(assembly_csi)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler),
          val("vamb"), path("vamb"), path("vamb.contig2bin.tsv"), path("vamb.note"),                         emit: result

    script:
    """
    run_vamb.sh \\
      --assembly "${assembly_fasta}" \\
      --bam "${assembly_bam}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries}
    """
}


process VAMB_GPU {
    tag "${sra}:${srr}:${assembler}"
    label 'gpu'
    publishDir { "${assemblerPublishDir(sra, srr, assembler)}/binning" },
      mode: 'copy',
      overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(assembly_fasta), path(assembly_bam), path(assembly_csi)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler),
          val("vamb"), path("vamb"), path("vamb.contig2bin.tsv"), path("vamb.note"),                         emit: result

    script:
    """
    run_vamb.sh \\
      --assembly "${assembly_fasta}" \\
      --bam "${assembly_bam}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries} \\
      --cuda \\
      --require-cuda
    """
}


process LORBIN {
    tag "${sra}:${srr}:${assembler}"
    label 'binning'
    publishDir { "${assemblerPublishDir(sra, srr, assembler)}/binning" },
      mode: 'copy',
      overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(assembly_fasta), path(assembly_bam), path(assembly_csi)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler),
          val("lorbin"), path("lorbin"), path("lorbin.contig2bin.tsv"), path("lorbin.note"),                  emit: result

    script:
    """
    run_lorbin.sh \\
      --assembly "${assembly_fasta}" \\
      --bam "${assembly_bam}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries}
    """
}


process LORBIN_GPU {
    tag "${sra}:${srr}:${assembler}"
    label 'gpu'
    publishDir { "${assemblerPublishDir(sra, srr, assembler)}/binning" },
      mode: 'copy',
      overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(assembly_fasta), path(assembly_bam), path(assembly_csi)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler),
          val("lorbin"), path("lorbin"), path("lorbin.contig2bin.tsv"), path("lorbin.note"),                  emit: result

    script:
    """
    run_lorbin.sh \\
      --assembly "${assembly_fasta}" \\
      --bam "${assembly_bam}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries} \\
      --require-cuda
    """
}


process SEMIBIN {
    tag "${sra}:${srr}:${assembler}"
    label 'binning'
    publishDir { "${assemblerPublishDir(sra, srr, assembler)}/binning" },
      mode: 'copy',
      overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(assembly_fasta), path(assembly_bam), path(assembly_csi)
    path uniprot_db

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler),
          val("semibin"), path("semibin"), path("semibin.contig2bin.tsv"), path("semibin.note"),             emit: result

    script:
    """
    run_semibin.sh \\
      --assembly "${assembly_fasta}" \\
      --bam "${assembly_bam}" \\
      --diamond-db "${uniprot_db}" \\
      --read-type "${read_type}" \\
      --environment "${params.semibin_environment}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries}
    """
}


process ROSELLA {
    tag "${sra}:${srr}:${assembler}"
    label 'binning'
    publishDir { "${assemblerPublishDir(sra, srr, assembler)}/binning" },
      mode: 'copy',
      overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path(assembly_fasta), path(assembly_bam), path(assembly_csi)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler),
          val("rosella"), path("rosella"), path("rosella.contig2bin.tsv"), path("rosella.note"),             emit: result

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
    tag "${sra}:${srr}:${assembler}"
    publishDir { "${assemblerPublishDir(sra, srr, assembler)}/binning" },
      mode: 'copy',
      overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler),
          path(assembly_fasta),
          path(contig2bin_maps)

    output:
    tuple val(sra), val(srr), val(read_type), val(assembler), path("dastool"),                                               emit: bins
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("dastool.note"), emit: note

    script:
    def mapArg = contig2bin_maps.collect { contig2bin_map -> contig2bin_map.toString() }.join(',')
    """
    run_dastool.sh \\
      --assembly "${assembly_fasta}" \\
      --bin-maps "${mapArg}" \\
      --cpus ${task.cpus} \\
      --attempt ${task.attempt} \\
      --max-retries ${params.max_retries}
    """
}


process BINETTE {
    tag "${sra}:${srr}:${assembler}"
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
}


process CREATE_EMPTY_SUMMARY {
  tag "${sra}:${srr}:${read_type}:${assembler}"

  input:
  tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), val(note)

  output:
  tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("empty_summary.csv"), val(note), emit: skipped_rows

  script:
  """
  echo 'rank,ncbi_taxa,n_contigs,output_ids_csv,output_fasta' > empty_summary.csv
  """
}


process CREATE_ASSEMBLER_SELECTION_NOTE {
  tag "${sra}:${srr}:${read_type}"

  input:
  tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), val(note)

  output:
  tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(read_type), val(assembler), path("FAIL.note"), emit: note

  script:
  """
  printf '%s\n' "${note}" > FAIL.note
  """
}


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


//-- Workflow ------------------------------------------------------------------
workflow PRE_SCREENING {
  take:
    sra_accessions_channel
    validated_taxa
    sandpiper_db_ch
    singlem_db_ch

  main:
    def doScreening = params.taxa != null

    // Start PRE_SCREENING after taxa validation
    def gated_sra_accessions_channel = doScreening
                        ? sra_accessions_channel.combine(validated_taxa).map { sra, vt -> sra }
                        : sra_accessions_channel

    // Extract metadata & filter SRR
    sra_metadata = DOWNLOAD_SRA_METADATA(gated_sra_accessions_channel)
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
                  def read_type = (row.read_type ?: row.assembler ?: '').trim()
                  [sra, srr, platform, model, strategy, read_type]
              }
              .filter { row_values -> row_values[1] }  // ensure SRR not empty
              .distinct()

    def singlem_reads_ch
    def sandpiper_note_ch
    def download_srr_note_ch
    def singlem_note_ch

    if (doScreening) {
      // SANDPIPER prescreening
      sandpiper = SANDPIPER(srr_ch, validated_taxa, sandpiper_db_ch)
      // decision: NEGATIVE / RUN_SINGLEM / PASS
      sandpiper_decision_ch = sandpiper.decision.map { sra, srr, platform, model, strategy, read_type, dec_path ->
        def decision = file(dec_path).text.trim()
        tuple(sra, srr, platform, model, strategy, read_type, decision)
      }

      srr_prescreened = sandpiper_decision_ch
        .filter { sra, srr, platform, model, strategy, read_type, decision ->
          decision == 'PASS' || decision == 'RUN_SINGLEM'
        }

      // Download SRR reads
      download_srr = DOWNLOAD_SRR(srr_prescreened)

      srr_reads = download_srr.reads.map { sra, srr, platform, model, strategy, detected_read_type, sandpiper_dec, reads, asm_txt ->
        def fixedReadType = file(asm_txt)?.text?.trim() ?: detected_read_type
        tuple(sra, srr, platform, model, strategy, fixedReadType, sandpiper_dec, reads)
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
      srr_prescreened = srr_ch.map { sra, srr, platform, model, strategy, read_type ->
        tuple(sra, srr, platform, model, strategy, read_type, 'PASS')
      }

      // Download SRR reads
      download_srr = DOWNLOAD_SRR(srr_prescreened)

      // Normalise read type, then drop sandpiper decision.
      singlem_reads_ch = download_srr.reads.map { sra, srr, platform, model, strategy, detected_read_type, sandpiper_dec, reads, asm_txt ->
        def fixedReadType = file(asm_txt)?.text?.trim() ?: detected_read_type
        tuple(sra, srr, platform, model, strategy, fixedReadType, reads)
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
    selected_reads = singlem_reads.flatMap { sra, srr, platform, model, strategy, read_type, reads ->
      assemblersForReadType(read_type).collect { assembler ->
        tuple(sra, srr, platform, model, strategy, read_type, assembler, reads)
      }
    }

    assembler_selection_errors = singlem_reads.flatMap { sra, srr, platform, model, strategy, read_type, reads ->
      if (assemblersForReadType(read_type)) {
        return []
      }
      def note = "no compatible assembler selected for read_type '${read_type}'"
      [tuple(sra, srr, platform, model, strategy, read_type, '', note)]
    }

    assembler_selection_note = CREATE_ASSEMBLER_SELECTION_NOTE(assembler_selection_errors)

    metaspades_reads_channel = selected_reads.filter { sra, srr, platform, model, strategy, read_type, assembler, reads ->
      read_type.equalsIgnoreCase('short') && assembler == 'metaspades'
    }
    unicycler_reads_channel = selected_reads.filter { sra, srr, platform, model, strategy, read_type, assembler, reads ->
      read_type.equalsIgnoreCase('short') && assembler == 'unicycler'
    }
    metaflye_nanopore_reads_channel = selected_reads.filter { sra, srr, platform, model, strategy, read_type, assembler, reads ->
      read_type.equalsIgnoreCase('nanopore') && assembler == 'metaflye'
    }
    metaflye_pacbio_reads_channel = selected_reads.filter { sra, srr, platform, model, strategy, read_type, assembler, reads ->
      read_type.equalsIgnoreCase('pacbio') && assembler == 'metaflye'
    }
    metaflye_hifi_reads_channel = selected_reads.filter { sra, srr, platform, model, strategy, read_type, assembler, reads ->
      read_type.equalsIgnoreCase('hifi') && assembler == 'metaflye'
    }
    myloasm_hifi_reads_channel = selected_reads.filter { sra, srr, platform, model, strategy, read_type, assembler, reads ->
      read_type.equalsIgnoreCase('hifi') && assembler == 'myloasm'
    }

    short_reads_channel = channel.empty()
                        .mix(metaspades_reads_channel)
                        .mix(unicycler_reads_channel)

    fastp_short_out = FASTP_SHORT(short_reads_channel)

    metaspades_preprocessed_reads_channel = fastp_short_out.reads.filter { sra, srr, platform, model, strategy, read_type, assembler, reads ->
      read_type.equalsIgnoreCase('short') && assembler == 'metaspades'
    }
    unicycler_preprocessed_reads_channel = fastp_short_out.reads.filter { sra, srr, platform, model, strategy, read_type, assembler, reads ->
      read_type.equalsIgnoreCase('short') && assembler == 'unicycler'
    }

    metaspades_out         = METASPADES(metaspades_preprocessed_reads_channel)
    unicycler_out          = UNICYCLER(unicycler_preprocessed_reads_channel)
    metaflye_nanopore_out  = METAFLYE_NANO(metaflye_nanopore_reads_channel)
    metaflye_pacbio_out    = METAFLYE_PACBIO(metaflye_pacbio_reads_channel)
    metaflye_hifi_out      = METAFLYE_HIFI(metaflye_hifi_reads_channel)
    myloasm_out            = MYLOASM(myloasm_hifi_reads_channel)

    // Step 6: DIAMOND
    assembly_fasta_channel = channel.empty()
                        .mix(metaspades_out.assembly_fasta)
                        .mix(unicycler_out.assembly_fasta)
                        .mix(metaflye_nanopore_out.assembly_fasta)
                        .mix(metaflye_pacbio_out.assembly_fasta)
                        .mix(metaflye_hifi_out.assembly_fasta)
                        .mix(myloasm_out.assembly_fasta)

    diamond = DIAMOND(assembly_fasta_channel, uniprot_db_ch)

    // Map reads back to assemblies for BlobTools and binning.
    mapping_reads_channel = channel.empty()
                  .mix(fastp_short_out.reads)
                  .mix(metaflye_nanopore_reads_channel)
                  .mix(metaflye_pacbio_reads_channel)
                  .mix(metaflye_hifi_reads_channel)
                  .mix(myloasm_hifi_reads_channel)

    assembly_fasta_for_mapping_by_sample_key = assembly_fasta_channel.map { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta ->
      tuple([sra, srr, read_type, assembler], [platform, model, strategy, assembly_fasta])
    }
    mapping_reads_by_sample_key = mapping_reads_channel.map { sra, srr, platform, model, strategy, read_type, assembler, reads ->
      tuple([sra, srr, read_type, assembler], [platform, model, strategy, reads])
    }
    mapping_input_channel = assembly_fasta_for_mapping_by_sample_key
      .join(mapping_reads_by_sample_key)
      .map { sample_join_key, assembly_payload, reads_payload ->
        def (sra, srr, read_type, assembler) = sample_join_key
        def (platform, model, strategy, assembly_fasta) = assembly_payload
        def reads = reads_payload[3]
        tuple(sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, reads)
      }

    map_to_assembly_out = MAP_TO_ASSEMBLY(mapping_input_channel)
    assembly_bam_channel = map_to_assembly_out.mapped

    // Step 7: BlobTools
    // Key every stream by sample, read type, and assembler tool.
    assembly_fasta_by_sample_key = assembly_fasta_channel.map  { sra, srr, platform, model, strategy, read_type, assembler, fasta -> tuple([sra, srr, read_type, assembler], [platform, model, strategy, fasta]) }
    diamond_blast_by_sample_key = diamond.blast.map { sra, srr, read_type, assembler, blast -> tuple([sra, srr, read_type, assembler], blast) }
    assembly_bam_by_sample_key   = assembly_bam_channel.map       { sra, srr, read_type, assembler, bam, csi -> tuple([sra, srr, read_type, assembler], [bam, csi]) }

    // Join (fasta * diamond) then * bam
    assembly_with_blast_by_sample_key     = assembly_fasta_by_sample_key.join(diamond_blast_by_sample_key)
    assembly_with_blast_and_bam_by_sample_key = assembly_with_blast_by_sample_key.join(assembly_bam_by_sample_key)

    // Unkey + call BlobTools
    blobtools_input_channel = assembly_with_blast_and_bam_by_sample_key.map { sample_join_key, assembly_payload, blast, alignment_files ->
      def (sra, srr, read_type, assembler) = sample_join_key
      def (platform, model, strategy, assembly) = assembly_payload
      def (bam, csi) = alignment_files
      tuple(sra, srr, platform, model, strategy, read_type, assembler, assembly, blast, bam, csi)
    }

    blobtools = BLOBTOOLS(blobtools_input_channel, taxdump_ch)

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
        .map { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, blobtable_csv ->
          // note field is empty string for successful runs
          tuple(sra, srr, platform, model, strategy, read_type, assembler, '')
        }
      no_taxa_summary = CREATE_EMPTY_SUMMARY(no_taxa_success_meta).skipped_rows

      taxa_summary_ch = no_taxa_summary.map { sra, srr, platform, model, strategy, read_type, assembler, summary_csv, note ->
        tuple(sra, srr, platform, model, strategy, read_type, assembler, summary_csv)
      }

      taxa_note_ch    = channel.empty()
    }

    assembly_notes_ch = channel.empty()
                          .mix(assembler_selection_note.note)
                          .mix(fastp_short_out.note)
                          .mix(metaspades_out.note)
                          .mix(unicycler_out.note)
                          .mix(metaflye_nanopore_out.note)
                          .mix(metaflye_pacbio_out.note)
                          .mix(metaflye_hifi_out.note)
                          .mix(myloasm_out.note)
                          .mix(map_to_assembly_out.note)

  emit:
    // For binning
    assembly_bam_all = assembly_bam_channel
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
    // blobtable_ch: (sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, blobtable)
    // assembly_bam_ch: (sra, srr, read_type, assembler, bam, csi)
    assembly_bam_by_sample_key = assembly_bam_ch.map { sra, srr, read_type, assembler, bam, csi ->
      tuple([sra, srr, read_type, assembler], [bam, csi])
    }
    binning_base_by = blobtable_ch.map { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, blobtable ->
      tuple([sra, srr, read_type, assembler], [platform, model, strategy, assembly_fasta])
    }
    binning_join = binning_base_by.join(assembly_bam_by_sample_key)

    // binning_input: (sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, bam, csi)
    binning_input = binning_join.map { sample_join_key, assembly_payload, bam_payload ->
      def (sra, srr, read_type, assembler) = sample_join_key
      def (platform, model, strategy, assembly_fasta) = assembly_payload
      def (bam, csi) = bam_payload
      tuple(sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, bam, csi)
    }

    def binning_options = validateBinningOptions()
    def selected_binners = binning_options.binners
    def selected_refiners = binning_options.refiners
    def use_gpu_binners = binning_options.gpu

    if (use_gpu_binners) {
      def cpu_only_binners = selected_binners.findAll { binner -> binner in ['metabat', 'rosella', 'semibin'] }
      if (cpu_only_binners) {
        log.warn "GPU mode requested; CPU-only binner(s) will remain on CPU: ${cpu_only_binners.join(', ')}"
      }
    }

    planned_binning_input = binning_input
      .map { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, bam, csi ->
        def sample_binner_csv = binnersForSample(binning_options, sra, srr, read_type).join(',')
        tuple(sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, bam, csi, sample_binner_csv)
      }

    binner_plan_by_sample = planned_binning_input
      .map { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, bam, csi, sample_binner_csv ->
        tuple([sra, srr, read_type, assembler], binnerCsvSize(sample_binner_csv))
      }

    metabat_input = planned_binning_input
      .filter { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, bam, csi, sample_binner_csv ->
        binnerCsvContains(sample_binner_csv, 'metabat')
      }
      .map { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, bam, csi, sample_binner_csv ->
        tuple(sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, bam, csi)
      }

    comebin_input = planned_binning_input
      .filter { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, bam, csi, sample_binner_csv ->
        binnerCsvContains(sample_binner_csv, 'comebin')
      }
      .map { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, bam, csi, sample_binner_csv ->
        tuple(sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, bam, csi)
      }

    vamb_input = planned_binning_input
      .filter { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, bam, csi, sample_binner_csv ->
        binnerCsvContains(sample_binner_csv, 'vamb')
      }
      .map { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, bam, csi, sample_binner_csv ->
        tuple(sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, bam, csi)
      }

    lorbin_input = planned_binning_input
      .filter { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, bam, csi, sample_binner_csv ->
        binnerCsvContains(sample_binner_csv, 'lorbin')
      }
      .map { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, bam, csi, sample_binner_csv ->
        tuple(sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, bam, csi)
      }

    semibin_input = planned_binning_input
      .filter { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, bam, csi, sample_binner_csv ->
        binnerCsvContains(sample_binner_csv, 'semibin')
      }
      .map { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, bam, csi, sample_binner_csv ->
        tuple(sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, bam, csi)
      }

    rosella_input = planned_binning_input
      .filter { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, bam, csi, sample_binner_csv ->
        binnerCsvContains(sample_binner_csv, 'rosella')
      }
      .map { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, bam, csi, sample_binner_csv ->
        tuple(sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, bam, csi)
      }

    metabat_results = channel.empty()
    comebin_results = channel.empty()
    vamb_results = channel.empty()
    lorbin_results = channel.empty()
    semibin_results = channel.empty()
    rosella_results = channel.empty()

    if ('metabat' in selected_binners) {
      metabat_results = METABAT(metabat_input).result
    }
    if ('comebin' in selected_binners) {
      if (use_gpu_binners) {
        comebin_results = COMEBIN_GPU(comebin_input).result
      }
      else {
        comebin_results = COMEBIN(comebin_input).result
      }
    }
    if ('vamb' in selected_binners) {
      if (use_gpu_binners) {
        vamb_results = VAMB_GPU(vamb_input).result
      }
      else {
        vamb_results = VAMB(vamb_input).result
      }
    }
    if ('lorbin' in selected_binners) {
      if (use_gpu_binners) {
        lorbin_results = LORBIN_GPU(lorbin_input).result
      }
      else {
        lorbin_results = LORBIN(lorbin_input).result
      }
    }
    if ('semibin' in selected_binners) {
      semibin_results = SEMIBIN(semibin_input, uniprot_db_ch).result
    }
    if ('rosella' in selected_binners) {
      rosella_results = ROSELLA(rosella_input).result
    }

    binner_results = channel.empty()
      .mix(metabat_results)
      .mix(comebin_results)
      .mix(vamb_results)
      .mix(lorbin_results)
      .mix(semibin_results)
      .mix(rosella_results)

    dastool_base_by = binning_input.map { sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta, bam, csi ->
      tuple([sra, srr, read_type, assembler], [platform, model, strategy, assembly_fasta])
    }

    binner_maps_by_sample = binner_results
      .map { sra, srr, platform, model, strategy, read_type, assembler, tool, bin_dir, contig2bin, note_path ->
        tuple([sra, srr, read_type, assembler], [tool, contig2bin])
      }
      .combine(binner_plan_by_sample, by: 0)
      .map { sample_join_key, binner_map_payload, expectedCount -> tuple(groupKey(sample_join_key, expectedCount as int), binner_map_payload) }
      .groupTuple()
      .map { grouped_sample_key, binner_map_entries -> tuple(grouped_sample_key.target, binner_map_entries) }

    dastool_join = dastool_base_by.join(binner_maps_by_sample)

    dastool_in = dastool_join.map { sample_join_key, assembly_payload, binner_map_entries ->
      def (sra, srr, read_type, assembler) = sample_join_key
      def (platform, model, strategy, assembly_fasta) = assembly_payload
      def contig2bin_maps = binner_map_entries.collect { binner_map_entry -> binner_map_entry[1] }
      tuple(
        sra, srr, platform, model, strategy, read_type, assembler, assembly_fasta,
        contig2bin_maps
      )
    }
    .filter { row -> row != null }

    dastool_note_entries = channel.empty()
    if ('dastool' in selected_refiners) {
      dastool_note_entries = DASTOOL(dastool_in).note.map { sra, srr, platform, model, strategy, read_type, assembler, note_path ->
        tuple(sra, srr, platform, model, strategy, read_type, assembler, 'dastool', note_path)
      }
    }

    binette_note_entries = channel.empty()
    if ('binette' in selected_refiners) {
      def checkm2_db_ch = channel.value(file(params.checkm2_db, checkIfExists: true))
      binette_note_entries = BINETTE(dastool_in, checkm2_db_ch).note.map { sra, srr, platform, model, strategy, read_type, assembler, note_path ->
        tuple(sra, srr, platform, model, strategy, read_type, assembler, 'binette', note_path)
      }
    }

    binner_note_entries = binner_results.map { sra, srr, platform, model, strategy, read_type, assembler, tool, bin_dir, contig2bin, note_path ->
      tuple(sra, srr, platform, model, strategy, read_type, assembler, tool, note_path)
    }

    binning_note_entries = channel.empty()
      .mix(binner_note_entries)
      .mix(dastool_note_entries)
      .mix(binette_note_entries)

  emit:
    note_entries = binning_note_entries
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

    // from BINNING (empty when binning disabled)
    binning_note_entries

    // constant
    outdir

  main:
    def doScreening = (params.taxa != null)
    def noAssembly  = params.noassembly?.toString()?.toBoolean() ?: false

    // EXTRACT_TAXA only exists when assembly is running
    def doExtraction = doScreening && !noAssembly

    // Binning only exists when assembly is running.
    def binningRequested = params.binning?.toString()?.toBoolean() ?: false
    def doBinning        = binningRequested && !noAssembly


    // 1) Convert skipped SRA metadata into per‑SRR rows with a textual note
    skipped_srr_rows_channel = sra_metadata_skipped
      .map { sra, csvfile -> file(csvfile) }
      .splitCsv(header: true, strip: true)
      .map { row ->
        def sra       = (row.accession           ?: '').trim()
        def srr       = (row.run_accession       ?: '').trim()
        def platform  = (row.instrument_platform ?: '').trim()
        def model     = (row.instrument_model    ?: '').trim()
        def strategy  = (row.library_strategy    ?: '').trim()
        def note      = "did not match the criteria: ${(row.skip_reason ?: '').trim()}"
        // Read type and assembler are empty for skipped rows.
        tuple(sra, srr, platform, model, strategy, '', '', note)
      }
      .filter { row_values -> row_values[1] } // keep only rows with srr


    // 2) Base "successful" samples: anything with a taxa_summary (summary.csv)
    //    Note starts empty here; will add taxa/binning notes later.
    successful_summary_rows_channel = taxa_summary.map { sra, srr, platform, model, strategy, read_type, assembler, summary_csv ->
      tuple(sra, srr, platform, model, strategy, read_type, assembler, summary_csv, '')
    }


    // 3) Attach taxa notes and classify soft vs fatal EXTRACT_TAXA outcomes
    def summary_rows_with_taxa_notes_channel = successful_summary_rows_channel
    def taxa_failure_rows_channel = channel.empty()

    if( doExtraction ) {
      // Exactly 2 entries per sample:
      // - [summary_csv, base_note] from succeeded_sra
      // - note_path from taxa_note
      def taxa_group_size = 2

      successful_summary_by_taxa_key = successful_summary_rows_channel.map { sra, srr, platform, model, strategy, read_type, assembler, summary_csv, base_note ->
        def sample_key_map = [sra: sra, srr: srr, platform: platform,
                              model: model, strategy: strategy, read_type: read_type, assembler: assembler]
        def grouped_sample_key = groupKey(sample_key_map, taxa_group_size)
        tuple(grouped_sample_key, [summary_csv, base_note])
      }

      taxa_note_by_sample_key = taxa_note.map { sra, srr, platform, model, strategy, read_type, assembler, note_path ->
        def sample_key_map = [sra: sra, srr: srr, platform: platform,
                              model: model, strategy: strategy, read_type: read_type, assembler: assembler]
        def grouped_sample_key = groupKey(sample_key_map, taxa_group_size)
        tuple(grouped_sample_key, note_path)
      }

      success_and_taxa_notes_channel = channel.empty()
        .mix(successful_summary_by_taxa_key)
        .mix(taxa_note_by_sample_key)
        .groupTuple()

      // Unpack to: (sra, srr, platform, model, strategy, read_type, assembler, summary_csv, base_note, taxa_text)
      summary_rows_with_taxa_text_channel = success_and_taxa_notes_channel
        .map { grouped_sample_key, grouped_payloads ->
          def sample_key = grouped_sample_key.target as Map
          def (sra, srr, platform, model, strategy, read_type, assembler) =
            [sample_key.sra, sample_key.srr, sample_key.platform, sample_key.model, sample_key.strategy, sample_key.read_type, sample_key.assembler]

          // [summary_csv, base_note]
          def summary_payload = grouped_payloads.find { payload -> payload instanceof List && payload.size() == 2 }
          if( !summary_payload ) {
            return null
          }

          def summary_csv = summary_payload[0]
          def base_note   = summary_payload[1] ?: ''

          // note_path for EXTRACT_TAXA
          def taxa_note_path = grouped_payloads.find { payload -> !(payload instanceof List) }
          def taxa_text = ''
          if( taxa_note_path ) {
            taxa_text = file(taxa_note_path as String).text.trim()
          }

          tuple(sra, srr, platform, model, strategy, read_type, assembler, summary_csv, base_note, taxa_text)
        }
        .filter { row -> row != null }

      // Split EXTRACT_TAXA outcomes into "success" vs "fatal"
      def softPattern = 'skipping extraction because taxa list contains only GTDB-style taxa'

      def taxa_branches = summary_rows_with_taxa_text_channel.branch { sra, srr, platform, model, strategy, read_type, assembler, summary_csv, base_note, taxa_text ->
        // Empty note or GTDB-only “soft skip” => success
        success: (!taxa_text || taxa_text.contains(softPattern))
        // Anything else => fatal
        fatal:   (taxa_text && !taxa_text.contains(softPattern))
      }

      taxa_success = taxa_branches.success
      taxa_fatal   = taxa_branches.fatal

      // Non-fatal EXTRACT_TAXA results remain as "successful" samples.
      // Append taxa_text (if any) to the note field.
      summary_rows_with_taxa_notes_channel = taxa_success.map { sra, srr, platform, model, strategy, read_type, assembler, summary_csv, base_note, taxa_text ->
        def final_note = base_note
        if( taxa_text ) {
          final_note = final_note ? "${final_note}; ${taxa_text}" : taxa_text
        }
        tuple(sra, srr, platform, model, strategy, read_type, assembler, summary_csv, final_note)
      }

      // Fatal EXTRACT_TAXA outcomes become "errors" that will go through
      // CREATE_EMPTY_SUMMARY, similar to other fatal notes.
      taxa_failure_rows_channel = taxa_fatal.map { sra, srr, platform, model, strategy, read_type, assembler, summary_csv, base_note, taxa_text ->
        def note_text = base_note
        if( taxa_text ) {
          note_text = note_text ? "${note_text}; ${taxa_text}" : taxa_text
        }
        tuple(sra, srr, platform, model, strategy, read_type, assembler, note_text)
      }
    }


    // 4) Aggregate binning notes into a single string per sample (if binning)
    def binning_notes_by_sample_channel = channel.empty()

    if( doBinning ) {
      grouped_binning_notes_channel = binning_note_entries
        .map { sra, srr, platform, model, strategy, read_type, assembler, tool, note_path ->
          def sample_key = [sra: sra, srr: srr, platform: platform, model: model,
                            strategy: strategy, read_type: read_type, assembler: assembler]
          tuple(sample_key, [tool, note_path])
        }
        .groupTuple()

      // binning_agg: (sra, srr, platform, model, strategy, read_type, assembler, "tool1: msg; tool2: msg; ...")
      binning_notes_by_sample_channel = grouped_binning_notes_channel.map { grouped_sample_key, note_entries ->
        def sample_key = grouped_sample_key as Map
        def (sra, srr, platform, model, strategy, read_type, assembler) =
          [sample_key.sra, sample_key.srr, sample_key.platform, sample_key.model, sample_key.strategy, sample_key.read_type, sample_key.assembler]

        def note_texts = note_entries.collect { note_entry ->
          def (tool, note_path) = note_entry
          def note_text = file(note_path as String).text.trim()
          note_text ? "${tool}: ${note_text}" : null
        }.findAll { token -> token }

        def joined_note_text = note_texts ? note_texts.join('; ') : ''
        tuple(sra, srr, platform, model, strategy, read_type, assembler, joined_note_text)
      }
    }


    // 5) Combine successful rows with aggregated binning notes (if binning)
    def final_success_rows_channel = summary_rows_with_taxa_notes_channel

    if( doBinning ) {
      def success_group_size = 2

      successful_summary_by_binning_key = summary_rows_with_taxa_notes_channel.map { sra, srr, platform, model, strategy, read_type, assembler, summary_csv, base_note ->
        def sample_key_map = [sra: sra, srr: srr, platform: platform,
                              model: model, strategy: strategy, read_type: read_type, assembler: assembler]
        def grouped_sample_key = groupKey(sample_key_map, success_group_size)
        tuple(grouped_sample_key, [summary_csv, base_note])
      }

      binning_note_by_success_key = binning_notes_by_sample_channel.map { sra, srr, platform, model, strategy, read_type, assembler, binning_note ->
        def sample_key_map = [sra: sra, srr: srr, platform: platform,
                              model: model, strategy: strategy, read_type: read_type, assembler: assembler]
        def grouped_sample_key = groupKey(sample_key_map, success_group_size)
        tuple(grouped_sample_key, binning_note)
      }

      success_and_binning_notes_channel = channel.empty()
        .mix(successful_summary_by_binning_key)
        .mix(binning_note_by_success_key)
        .groupTuple()

      final_success_rows_channel = success_and_binning_notes_channel
        .map { grouped_sample_key, grouped_payloads ->
          def sample_key = grouped_sample_key.target as Map
          def (sra, srr, platform, model, strategy, read_type, assembler) =
            [sample_key.sra, sample_key.srr, sample_key.platform, sample_key.model, sample_key.strategy, sample_key.read_type, sample_key.assembler]

          def summary_payload = grouped_payloads.find { payload -> payload instanceof List && payload.size() == 2 }
          if( !summary_payload ) {
            return null
          }

          def summary_csv = summary_payload[0]
          def base_note   = summary_payload[1] ?: ''

          def binning_note_text = ''
          def extra_note_text = grouped_payloads.find { payload -> !(payload instanceof List) } as String
          if( extra_note_text ) {
            binning_note_text = extra_note_text.trim()
          }

          def final_note = base_note
          if( binning_note_text ) {
            final_note = final_note ? "${final_note}; ${binning_note_text}" : binning_note_text
          }

          tuple(sra, srr, platform, model, strategy, read_type, assembler, summary_csv, final_note)
        }
        .filter { row -> row != null }
    }


    // 6) Collect all fatal errors (including fatal EXTRACT_TAXA) and turn them
    //    into empty per-sample summaries.
    failure_rows_channel = channel.empty()
      .mix(sra_metadata_note)
      .mix(sandpiper_note)
      .mix(download_srr_note)
      .mix(singlem_note)
      .mix(assembly_notes)
      .mix(diamond_note)
      .mix(blobtools_note)
      .map { sra, srr, platform, model, strategy, read_type, assembler, note_path ->
        def note = file(note_path).text.trim()
        tuple(sra, srr, platform, model, strategy, read_type, assembler, note)
      }
      // Fatal EXTRACT_TAXA errors (inc. "run failed") are added here
      .mix(taxa_failure_rows_channel)
      // Plus skipped SRR rows
      .mix(skipped_srr_rows_channel)

    // failed_sra: (sra, srr, platform, model, strategy, read_type, assembler, empty_summary.csv, note)
    failed_summary_rows_channel = CREATE_EMPTY_SUMMARY(failure_rows_channel)


    // 7) Combine succeeded and failed, and append rows to summary.tsv
    summary_rows_channel = channel.empty()
      .mix(final_success_rows_channel)
      .mix(failed_summary_rows_channel)

    summary_result = APPEND_SUMMARY(summary_rows_channel, outdir)

  emit:
    global_summary = summary_result.global_summary
}


workflow {
    main:
    if (params.help) {
      helpMessage()
      exit 0
    }

    if (params.gpu?.toString()?.toBoolean()) {
      error "GPU mode is planned but not implemented yet for the current binning selection workflow"
    }

    def sraMode   = params.sra != null
    def fastqMode = params.fastq_tsv != null

    if (!sraMode && !fastqMode) {
      log.error "Error: either --sra or --fastq_tsv parameter must be provided"
      missingParametersError()
    }

    def doScreening = params.taxa != null
    def noAssembly = params.noassembly?.toString()?.toBoolean() ?: false
    def doAssembly = !noAssembly
    if (doAssembly) {
      selectedAssemblerTokens()
    }

    // If noassembly, binning makes no sense.
    def binningRequested = params.binning?.toString()?.toBoolean() ?: false
    def doBinning        = doAssembly && binningRequested
    if (!doAssembly && binningRequested) {
      log.warn "Warning: --binning is ignored because --noassembly was set"
    }

    // Only require assembly params when assembly is enabled
    if (doAssembly) {
      if (!params.taxdump || !params.uniprot_db) {
        log.error "Error: Missing --taxdump or --uniprot_db"
        missingParametersError()
      }
    }

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
    else if (!doAssembly) {
      log.warn "Warning: --noassembly was set but --taxa was not provided; SINGLEM will not run in this mode."
    }


    // Common channels
    def outdir = file(params.outdir ?: './output').toAbsolutePath().toString()

    // taxdump is needed for screening and/or assembly; keep it optional when not needed
    def taxdump_ch = params.taxdump ? channel.value(file(params.taxdump)) : channel.empty()

    // uniprot is only needed for assembly (DIAMOND/SEMIBIN)
    def uniprot_db_ch = (doAssembly && params.uniprot_db) ? channel.value(file(params.uniprot_db)) : channel.empty()


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
      sra_accessions_channel = channel.fromPath(params.sra, checkIfExists: true)
                      .splitCsv(header: true, strip: true)
                      .map { row -> row.sra.trim() }
                      .filter { row -> row }
                      .distinct()

      pre_screening_out = PRE_SCREENING(sra_accessions_channel, validated_taxa_ch, sandpiper_db_ch, singlem_db_ch)

      sra_singlem_reads_ch    = pre_screening_out.singlem_reads
      sra_metadata_skipped_ch = pre_screening_out.sra_metadata_skipped
      sra_metadata_note       = pre_screening_out.sra_metadata_note
      sra_sandpiper_note      = pre_screening_out.sandpiper_note
      sra_download_srr_note   = pre_screening_out.download_srr_note
      sra_singlem_note        = pre_screening_out.singlem_note
    }

    // FASTQ mode
    def fastq_singlem_reads_ch = channel.empty()
    def fastq_singlem_note     = channel.empty()
    if (fastqMode) {
      fastq_samplesheet_channel = channel.fromPath(params.fastq_tsv, checkIfExists: true)
                        .splitCsv(header: true, sep: '\t', strip: true)
                        .map { row ->
                          def sample    = (row.sample ?: '').trim()
                          def read_type = (row.read_type ?: '').trim()
                          def reads_raw = (row.reads ?: '').trim()

                          if (!sample || !read_type || !reads_raw) {
                            log.warn "Skipping FASTQ TSV row with missing fields: ${row}"
                            return null
                          }

                          def read_files = reads_raw.split(/\s*,\s*/).findAll { token -> token }.collect { read_path -> file(read_path) }

                          if (!read_files) {
                            log.warn "No valid FASTQ paths for sample ${sample}; skipping"
                            return null
                          }

                          def sra       = sample
                          def srr       = sample
                          def platform  = "UNKNOWN"
                          def model     = read_type
                          def strategy  = "UNKNOWN"

                          tuple(sra, srr, platform, model, strategy, read_type, read_files)
                        }
                        .filter { row -> row != null }

      if (doScreening) {
        fastq_for_singlem_ch = fastq_samplesheet_channel.map { sra, srr, platform, model, strategy, read_type, reads ->
          tuple(sra, srr, platform, model, strategy, read_type, "RUN_SINGLEM", reads)
        }

        singlem_fastq = SINGLEM(fastq_for_singlem_ch, validated_taxa_ch, singlem_db_ch)
        fastq_singlem_reads_ch = singlem_fastq.reads
        fastq_singlem_note     = singlem_fastq.note
      }
      else {
        fastq_singlem_reads_ch = fastq_samplesheet_channel
        fastq_singlem_note     = channel.empty()
      }
    }

    // Merge SRA and FASTQ reads for assembly
    def screened_reads_channel = channel.empty()
                                  .mix(sra_singlem_reads_ch)
                                  .mix(fastq_singlem_reads_ch)
    def prescreening_note_channel = channel.empty()
                                  .mix(sra_singlem_note)
                                  .mix(fastq_singlem_note)


  // Empty notes
  def assembly_notes_ch = channel.empty()
  def diamond_note_ch   = channel.empty()
  def blobtools_note_ch = channel.empty()
  def taxa_note_ch      = channel.empty()
  def taxa_summary_ch   = channel.empty()
  def binning_note_entries_ch = channel.empty()

  if (noAssembly) {
    log.info "--noassembly set: skipping ASSEMBLY/BINNING; generating screening-only summary.tsv"

    // Prepare empty summary for succeeded samples
    def prescreen_success_meta = screened_reads_channel
      .map { sra, srr, platform, model, strategy, read_type, reads ->
        tuple(sra, srr, platform, model, strategy, read_type, '', '')
      }
      .distinct()
    def prescreen_empty = CREATE_EMPTY_SUMMARY(prescreen_success_meta).skipped_rows

    // Convert to the shape expected by SUMMARY's taxa_summary input:
    // (sra, srr, platform, model, strategy, read_type, assembler, summary_csv)
    taxa_summary_ch = prescreen_empty.map { sra, srr, platform, model, strategy, read_type, assembler, summary_csv, note ->
      tuple(sra, srr, platform, model, strategy, read_type, assembler, summary_csv)
    }
  }
  else {
    def assembly_out = ASSEMBLY(screened_reads_channel, validated_taxa_ch, uniprot_db_ch, taxdump_ch)

    assembly_notes_ch = assembly_out.assembly_notes
    diamond_note_ch   = assembly_out.diamond_note
    blobtools_note_ch = assembly_out.blobtools_note
    taxa_note_ch      = assembly_out.taxa_note
    taxa_summary_ch   = assembly_out.taxa_summary

    if (doBinning) {
      def binning_out = BINNING(assembly_out.blobtable, assembly_out.assembly_bam_all, uniprot_db_ch)
      binning_note_entries_ch = binning_out.note_entries
    }
  }


  SUMMARY(
    // PRE_SCREENING: summary-related outputs
    sra_metadata_skipped_ch,
    sra_metadata_note,
    sra_sandpiper_note,
    sra_download_srr_note,
    prescreening_note_channel,

    // ASSEMBLY: summary-related outputs (empty in --noassembly)
    assembly_notes_ch,
    diamond_note_ch,
    blobtools_note_ch,
    taxa_note_ch,
    taxa_summary_ch,

    // BINNING: raw note entries (empty in --noassembly or if binning off)
    binning_note_entries_ch,

    // constant
    outdir
  )

  workflow.onComplete = {
    def outdirPath  = file(params.outdir ?: './output').toAbsolutePath()
    def summaryFile = outdirPath.resolve('summary.tsv')
    def traceFile   = file("${workflow.launchDir}/execution-reports/trace.tsv").toAbsolutePath()
    def scriptFile  = file("${workflow.projectDir}/bin/annotate_summary_from_trace.py").toAbsolutePath()

    log.info "onComplete: summary.tsv -> ${summaryFile}"
    log.info "onComplete: trace.tsv   -> ${traceFile}"
    log.info "onComplete: annotator   -> ${scriptFile}"

    if( !summaryFile.exists() ) {
      log.warn "onComplete: ${summaryFile} not found; skipping scheduler annotation"
    }
    else if( !traceFile.exists() ) {
      log.warn "onComplete: ${traceFile} not found; skipping scheduler annotation"
    }
    else {
      def cmd = [
        'python3',
        scriptFile.toString(),
        summaryFile.toString(),
        traceFile.toString()
      ]

      log.info "onComplete: running ${cmd.join(' ')}"

      def proc = new java.lang.ProcessBuilder(cmd)
        .directory( workflow.launchDir.toFile() )
        .redirectError( java.lang.ProcessBuilder.Redirect.INHERIT )
        .redirectOutput( java.lang.ProcessBuilder.Redirect.INHERIT )
        .start()

      def rc = proc.waitFor()
      if( rc != 0 ) {
        log.warn "onComplete: annotator script exited with code ${rc}"
      }
      else {
        log.info "onComplete: summary.tsv successfully annotated with scheduler error information"
      }
    }
  }
}
