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
    --taxa          Taxa for extraction (e.g., phylum, genus)
    --taxdump       Path to taxdump json database
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
    error "Please provide all required parameters: --sra, --taxa, --taxdump, and --uniprot_db"
}


//-- Processes -----------------------------------------------------------------

process VALIDATE_TAXA {
    input:
    path(taxa_file)
    path(taxdump)

    output:
    path("validated_taxa.csv"), emit: valid_taxa

    script:
    """
    validate_taxa.py --taxa ${taxa_file} --taxdump ${taxdump} \\
      && cp -v ${taxa_file} validated_taxa.csv
    """
}


process DOWNLOAD_SRA_METADATA {
    tag "${sra}"
    publishDir "${params.outdir}/metadata/", mode: 'copy', overwrite: true

    input:
    val sra
    path(valid_taxa)

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


process DOWNLOAD_SRR {
    tag { "${sra}:${srr}" }

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("*.f*q*"),    optional:true, emit: reads
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("FAIL.note"), optional:true, emit: note

    script:
    """
    # Download sequence data
    if ! iseq -g -t ${task.cpus} -p 8 -i "${srr}"; then
      if [[ ${task.attempt} -lt ${params.max_retries} ]]; then
        exit 1
      fi
      echo "Download raw data failed" > FAIL.note
      rm -f *.f*q* "${srr}"  # remove sra and fastq files that didn't download properly
      exit 0
    fi
    """
}


process METASPADES {
    tag { "${sra}:${srr}" }
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
    R1=\$(ls *_1.fastq.gz || ls *_R1*.fastq.gz || true)
    R2=\$(ls *_2.fastq.gz || ls *_R2*.fastq.gz || true)
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
    tag { "${sra}:${srr}" }
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
    tag { "${sra}:${srr}" }
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
    tag { "${sra}:${srr}" }
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
    tag { "${sra}:${srr}" }
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
    tag { "${sra}:${srr}" }
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
    tag { "${sra}:${srr}" }
    publishDir "${params.outdir}/${sra}/${srr}/", mode: 'copy', overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path(assembly_fasta), path(blobtable)
    path(taxa)
    path(taxdump)

    output:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path("summary.csv"), emit: summary
    tuple val(sra), val(srr), path("*.ids.csv"), optional:true,emit: extracted_ids
    tuple val(sra), val(srr), path("*.fasta"), optional:true, emit: extracted_fasta

    script:
    """
    extract_records.py --blobtable "${blobtable}" --fasta "${assembly_fasta}" \\
      --taxa "${taxa}" --taxdump "${taxdump}"
    """
}


process APPEND_SUMMARY {
    tag { "${sra}:${srr}" }

    input:
    tuple val(sra), val(srr), val(platform), val(model), val(strategy), val(assembler), path(summary_csv)
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
    LINE="${sra}\t${srr}\t${platform}\t${model}\t${strategy}\t${assembler}\t\${COUNTS}\t"
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

  if (!params.sra || !params.uniprot_db || !params.taxa || !params.taxdump) {
    missingParametersError()
    exit 1
  }

  // Channel setup
  outdir = file(params.outdir).toAbsolutePath().toString()
  sra_ch = Channel.fromPath(params.sra, checkIfExists: true)
                  .splitCsv(header: true, strip: true)
                  .map { row -> row.sra.trim() }
                  .filter { it }
                  .distinct()
  taxa_ch    = Channel.value( file(params.taxa) )
  uniprot_db_ch = Channel.value( file(params.uniprot_db) )
  taxdump_ch    = Channel.value( file(params.taxdump) )

  // Validate taxa
  validated_taxa = VALIDATE_TAXA(taxa_ch, taxdump_ch)

  // Step 1: extract metadata & filter SRR
  sra_metadata = DOWNLOAD_SRA_METADATA(sra_ch, validated_taxa)
  filtered_srr = sra_metadata.filtered_sra
  skipped_srr  = sra_metadata.skipped_sra
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

  // Step 2: download SRR reads
  srr_reads = DOWNLOAD_SRR(srr_ch)
  srr_reads.reads.view { acc, srr, platform, model, strategy, asm, reads ->
    "${acc}\t${srr}\t\t${platform}\t${model}\t${strategy}\t${asm}\t${reads}"
  }

  // Step 3: assemble reads
  short_ch       = srr_reads.filter { sra, srr, platform, model, strategy, assembler, reads -> assembler.equalsIgnoreCase('short') }
  long_nano_ch   = srr_reads.filter { sra, srr, platform, model, strategy, assembler, reads -> assembler.equalsIgnoreCase('long_nano') }
  long_pacbio_ch = srr_reads.filter { sra, srr, platform, model, strategy, assembler, reads -> assembler.equalsIgnoreCase('long_pacbio') }
  long_hifi_ch   = srr_reads.filter { sra, srr, platform, model, strategy, assembler, reads -> assembler.equalsIgnoreCase('long_hifi') }

  spades_asm    = METASPADES(short_ch)
  flyenano_asm  = METAFLYE_NANO(long_nano_ch)
  flyepacbio_asm= METAFLYE_PACBIO(long_pacbio_ch)
  hifimeta_asm  = MYLOASM(long_hifi_ch)

  // Step 4: run DIAMOND
  asm_fasta_ch = Channel.empty()
                        .mix(spades_asm.assembly_fasta)
                        .mix(flyenano_asm.assembly_fasta)
                        .mix(flyepacbio_asm.assembly_fasta)
                        .mix(hifimeta_asm.assembly_fasta)


  diamond_ch = DIAMOND(asm_fasta_ch, uniprot_db_ch)
  // diamond_ch.view { sra, srr, asm, blast ->
  //   "DIAMOND: ${sra}\t${srr}\t${asm}\t${blast}"
  // }

  // Step 5. run BlobTools
  // Merge all BAM+Bai streams
  bam_ch = Channel.empty()
                  .mix(spades_asm.assembly_bam)
                  .mix(flyenano_asm.assembly_bam)
                  .mix(flyepacbio_asm.assembly_bam)
                  .mix(hifimeta_asm.assembly_bam)

  // Key every stream by (sra,srr,assembler)
  fasta_by  = asm_fasta_ch.map   { sra, srr, platform, model, strategy, assembler, fasta -> tuple([sra,srr], [platform,model,strategy,assembler,fasta]) }
  blast_by  = diamond_ch.map     { sra, srr, hits  -> tuple([sra,srr], hits) }
  bam_by    = bam_ch.map         { sra, srr, bam, bai -> tuple([sra,srr], [bam,bai]) }

  // Join (fasta * diamond) then * bam
  fasta_blast = fasta_by.join(blast_by)
  fasta_blast_bam = fasta_blast.join(bam_by)
  // -> ( [sra,srr,asm], fasta, hits, [bam,bai] )

  // Unkey + call BlobTools
  blobtools_in = fasta_blast_bam.map { key, fasta, hits, pair ->
    def (sra, srr) = key
    def (platform, model, strategy, assembler, assembly) = fasta
    def (bam, bai) = pair
    tuple(sra, srr, platform, model, strategy, assembler, assembly, hits, bam, bai)
  }

  blobtools_result = BLOBTOOLS(blobtools_in, taxdump_ch)
  // blobtools_result.view { sra, srr, asm, tbl -> "BLOBTOOLS:\t${sra}\t${srr}\t${asm}\t${tbl}" }

  // Step 6: extract taxa
  extracted_taxa = EXTRACT_TAXA(blobtools_result.blobtable, validated_taxa, taxdump_ch)

  // Step 7: append to global summary
  global_summary = APPEND_SUMMARY(extracted_taxa.summary, outdir)
}
