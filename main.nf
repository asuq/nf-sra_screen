#!/usr/bin/env nextflow
nextflow.enable.dsl=2


//-- Help Message ---------------------------------------------------------------

def helpMessage() {
  log.info """
  ===============================
  nf-sra_screen
  Nextflow pipeline for screening SRA genomes
  Version: 0.1.0
  Author : Akito Shima (ASUQ)
  Email: asuq.4096@gmail.com
  ===============================
  Usage: nextflow run main.nf [parameters]

  Required parameters:
    --sra           Path to sra.csv (header: sra)
    --nt_db         Path to BLAST nt database
    --blobtools_db  Path to BlobTools database

  Optional parameters:
    --help          Show this help message
    --outdir        Output directory (default: ./output)
  """.stripIndent()
}


def missingParametersError() {
    log.error "Missing input parameters"
    helpMessage()
    error "Please provide all required parameters: --sra, --nt_db, and --blobtools_db"
}


//-- Processes -----------------------------------------------------------------

process DOWNLOAD_SRA_METADATA {
    tag "${sra}"
    publishDir "${params.outdir}/metadata/${sra}", mode: 'copy', overwrite: true

    input:
    val sra

    output:
    tuple val(sra), path("${sra}.filtered.csv"), emit: filtered_sra

    script:
    """
    # Retrieve metadata
    iseq -m -i "${sra}"

    # Extract SRR to screen
    filter_sra.sh "${sra}.metadata.tsv"
    """
}

process DOWNLOAD_SRR {
    tag { "${sra}:${srr}" }
    //TODO: remove publishDir at the end to save storage
    publishDir "${params.outdir}/${sra}/${srr}", mode: 'copy', overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(strategy), val(model), val(assembler)

    output:
    tuple val(sra), val(srr), val(platform), val(strategy), val(model), val(assembler), path("*.f*q*"), emit: reads

    script:
    """
    # Download sequence data
    iseq -g -t ${task.cpus} -p 8 -i "${srr}"
    """
}


process METASPADES {
    tag { "${sra}:${srr}" }
    publishDir { "${params.outdir}/${sra}/${srr}/" }, mode: 'copy', overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(strategy), val(model), val(assembler), path(reads)

    output:
    tuple val(sra), val(srr), val(assembler), path("assembly.fasta"), emit: assembly_fasta
    tuple val(sra), val(srr), val(assembler), path("assembly.gfa"), emit: assembly_graph
    tuple val(sra), val(srr), val(assembler), path("assembly.log"), emit: assembly_log
    tuple val(sra), val(srr), val(assembler), path("fastp.html"), emit: fastp_html

    script:
    """
    R1=\$(ls *_1.fastq.gz || ls *_R1*.fastq.gz || true)
    R2=\$(ls *_2.fastq.gz || ls *_R2*.fastq.gz || true)
    if [[ -z "\$R1" || -z "\$R2" ]]; then
      echo "ERROR: Paired-end reads not found for ${srr}. Got R1='\$R1' R2='\$R2'" >&2
      exit 1
    fi

    # Run fastp
    fastp --in1 "\${R1}" --in2 "\${R2}" --out1 "${srr}_fastp_R1.fastq.gz" --out2 "${srr}_fastp_R2.fastq.gz" --thread ${task.cpus} \\
      --length_required 50 --detect_adapter_for_pe --html fastp.html

    # Run SPAdes
    spades.py -1 "${srr}_fastp_R1.fastq.gz" -2 "${srr}_fastp_R2.fastq.gz" \\
      -k 21,33,55,77,99,119,127 --meta --threads ${task.cpus} --memory ${task.memory.toGiga()}

    # Rename outputs
    mv -v scaffolds.fasta assembly.fasta
    mv -v assembly_graph_with_scaffolds.gfa assembly.gfa
    mv -v spades.log assembly.log
    touch assembly_info.txt
    """
}


process METAFLYE_NANO {
    tag { "${sra}:${srr}" }
    publishDir { "${params.outdir}/${sra}/${srr}/" }, mode: 'copy', overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(strategy), val(model), val(assembler), path(reads)

    output:
    tuple val(sra), val(srr), val(assembler), path("assembly.fasta"), emit: assembly_fasta
    tuple val(sra), val(srr), val(assembler), path("assembly.gfa"), emit: assembly_graph
    tuple val(sra), val(srr), val(assembler), path("assembly.log"), emit: assembly_log
    tuple val(sra), val(srr), val(assembler), path("assembly_info.txt"), emit: assembly_info

    script:
    """
    # Run metaFlye
    flye --nano-raw ${reads} --threads ${task.cpus} --scaffold --outdir '.' --meta

    # Rename outputs
    mv -v assembly_graph.gfa assembly.gfa
    mv -v flye.log assembly.log
    """
}


process METAFLYE_PACBIO {
    tag { "${sra}:${srr}" }
    publishDir { "${params.outdir}/${sra}/${srr}/" }, mode: 'copy', overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(strategy), val(model), val(assembler), path(reads)

    output:
    tuple val(sra), val(srr), val(assembler), path("assembly.fasta"), emit: assembly_fasta
    tuple val(sra), val(srr), val(assembler), path("assembly.gfa"), emit: assembly_graph
    tuple val(sra), val(srr), val(assembler), path("assembly.log"), emit: assembly_log
    tuple val(sra), val(srr), val(assembler), path("assembly_info.txt"), emit: assembly_info

    script:
    """
    # Run metaFlye
    flye --pacbio-raw ${reads} --threads ${task.cpus} --scaffold --outdir '.' --meta

    # Rename outputs
    mv -v assembly_graph.gfa assembly.gfa
    mv -v flye.log assembly.log
    """
}


process HIFIASM_META {
    tag { "${sra}:${srr}" }
    publishDir { "${params.outdir}/${sra}/${srr}/" }, mode: 'copy', overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(strategy), val(model), val(assembler), path(reads)

    output:
    tuple val(sra), val(srr), val(assembler), path("assembly.fasta"), emit: assembly_fasta
    tuple val(sra), val(srr), val(assembler), path("assembly.gfa"), emit: assembly_graph
    tuple val(sra), val(srr), val(assembler), path("assembly.log"), emit: assembly_log
    tuple val(sra), val(srr), val(assembler), path("assembly.bam"), emit: assembly_bam

    script:
    """
    # Run Hifimeta
    hifiasm_meta -t ${task.cpus} -o assembly ${reads} 2> assembly.log

    # Convert gfa file into fasta file
    awk '/^S/{print ">"\$2;print \$3}' "assembly.p_ctg.gfa" > "assembly.fasta"

    # Run minimap2
    minimap2 -ax map-hifi -t ${task.cpus} assembly.fasta ${reads} \\
      | samtools sort --output-fmt BAM -@ ${task.cpus} -o assembly.bam

    samtools index -b -@ ${task.cpus} assembly.bam

    # Rename outputs
    mv -v assembly.p_ctg.gfa assembly.gfa
    """
}


//-- Workflow ------------------------------------------------------------------
workflow {
  // Parameter parsing
  if (params.help) {
    helpMessage()
    exit 0
  }

  //TODO: Add check for nt_db and blobtools_db
  if (!params.sra) {
    missingParametersError()
    exit 1
  }

  // Channel setup
  sra_ch = Channel.fromPath(params.sra, checkIfExists: true)
                  .splitCsv(header: true, strip: true)
                  .map { row -> row.sra.trim() }
                  .filter { it }
                  .distinct()

  // Step 1: extract metadata & filter SRR
  filtered_srr = DOWNLOAD_SRA_METADATA(sra_ch)

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
                [sra, srr, platform, strategy, model, assembler]
            }
            .filter { it[1] }  // Ensure run_accession (SRR) is not empty
            .distinct()

// srr_ch.view { acc, srr, platform, strategy, model, asm ->
//       "DEBUG: ${acc}\t${srr}\t${platform}\t${strategy}\t${model}\t${asm}"
//   }

  // Step 2: download SRR reads
  srr_reads = DOWNLOAD_SRR(srr_ch)
  // srr_reads.view { acc, srr, platform, strategy, model, asm, reads ->
  //   "${acc}\t${srr}\t\t${platform}\t${strategy}\t${model}\t${asm}\t${reads}"
  // }

  // Step 3: assemble reads
  assembly = srr_reads.branch(
  short       : { sra, srr, platform, strategy, model, assembler, reads -> assembler.equalsIgnoreCase('short') },
  long_nano   : { sra, srr, platform, strategy, model, assembler, reads -> assembler.equalsIgnoreCase('long_nano') },
  long_pacbio : { sra, srr, platform, strategy, model, assembler, reads -> assembler.equalsIgnoreCase('long_pacbio') },
  long_hifi   : { sra, srr, platform, strategy, model, assembler, reads -> assembler.equalsIgnoreCase('long_hifi') }
  )

  spades_asm = METASPADES(assembly.short)
  flyenano_asm = METAFLYE_NANO(assembly.long_nano)
  flyepacbio_asm = METAFLYE_PACBIO(assembly.long_pacbio)
  hifimeta_asm = HIFIASM_META(assembly.long_hifi)
}
