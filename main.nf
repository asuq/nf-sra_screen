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
    publishDir "${params.outdir}/${sra}/${srr}/", mode: 'copy', overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(strategy), val(model), val(assembler), path(reads)

    output:
    tuple val(sra), val(srr), val(assembler), val(strategy), val(model), path("assembly.fasta"), emit: assembly_fasta
    tuple val(sra), val(srr), path("assembly.gfa"), emit: assembly_graph
    tuple val(sra), val(srr), path("assembly.log"), emit: assembly_log
    tuple val(sra), val(srr), path("assembly.bam"), path("assembly.bam.csi"), emit: assembly_bam
    tuple val(sra), val(srr), path("fastp.html"), emit: fastp_html

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
    spades.py -1 "${srr}_fastp_R1.fastq.gz" -2 "${srr}_fastp_R2.fastq.gz" -o '.' \\
      -k 21,33,55,77,99,119,127 --meta --threads ${task.cpus} --memory ${task.memory.toGiga()}

    # Rename outputs
    if [ -f scaffolds.fasta ]; then
      mv -v scaffolds.fasta assembly.fasta ;
    else
      mv -v contigs.fasta assembly.fasta ;
    fi
    mv -v assembly_graph_with_scaffolds.gfa assembly.gfa
    mv -v spades.log assembly.log

    # Run bowtie2
    bowtie2-build -f -q --threads ${task.cpus} assembly.fasta assembly_index
    bowtie2 -q --reorder --threads ${task.cpus} --time --met-stderr --met 10 \\
      -x assembly_index -1 "${srr}_fastp_R1.fastq.gz" -2 "${srr}_fastp_R2.fastq.gz" \\
      | samtools sort --output-fmt BAM -@ ${task.cpus} -o assembly.bam
    samtools index -c -o assembly.bam.csi -@ ${task.cpus} assembly.bam
    """
}


process METAFLYE_NANO {
    tag { "${sra}:${srr}" }
    publishDir "${params.outdir}/${sra}/${srr}/", mode: 'copy', overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(strategy), val(model), val(assembler), path(reads)

    output:
    tuple val(sra), val(srr), val(assembler), val(strategy), val(model), path("assembly.fasta"), emit: assembly_fasta
    tuple val(sra), val(srr), path("assembly.gfa"), emit: assembly_graph
    tuple val(sra), val(srr), path("assembly.log"), emit: assembly_log
    tuple val(sra), val(srr), path("assembly.bam"), path("assembly.bam.csi"), emit: assembly_bam

    script:
    """
    # Run metaFlye
    flye --nano-raw ${reads} --threads ${task.cpus} --scaffold --out-dir '.' --meta

    # Run minimap2
    minimap2 -ax map-ont -t ${task.cpus} assembly.fasta ${reads} \\
      | samtools sort --output-fmt BAM -@ ${task.cpus} -o assembly.bam
    samtools index -c -o assembly.bam.csi -@ ${task.cpus} assembly.bam

    # Rename outputs
    mv -v assembly_graph.gfa assembly.gfa
    mv -v flye.log assembly.log
    """
}


process METAFLYE_PACBIO {
    tag { "${sra}:${srr}" }
    publishDir "${params.outdir}/${sra}/${srr}/", mode: 'copy', overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(strategy), val(model), val(assembler), path(reads)

    output:
    tuple val(sra), val(srr), val(assembler), val(strategy), val(model), path("assembly.fasta"), emit: assembly_fasta
    tuple val(sra), val(srr), path("assembly.gfa"), emit: assembly_graph
    tuple val(sra), val(srr), path("assembly.log"), emit: assembly_log
    tuple val(sra), val(srr), path("assembly.bam"), path("assembly.bam.csi"), emit: assembly_bam

    script:
    """
    # Run metaFlye
    flye --pacbio-raw ${reads} --threads ${task.cpus} --scaffold --out-dir '.' --meta

    # Run minimap2
    minimap2 -ax map-pb -t ${task.cpus} assembly.fasta ${reads} \\
      | samtools sort --output-fmt BAM -@ ${task.cpus} -o assembly.bam
    samtools index -c -o assembly.bam.csi -@ ${task.cpus} assembly.bam

    # Rename outputs
    mv -v assembly_graph.gfa assembly.gfa
    mv -v flye.log assembly.log
    """
}

process MYLOASM {
    tag { "${sra}:${srr}" }
    publishDir "${params.outdir}/${sra}/${srr}/", mode: 'copy', overwrite: true

    input:
    tuple val(sra), val(srr), val(platform), val(strategy), val(model), val(assembler), path(reads)

    output:
    tuple val(sra), val(srr), val(assembler), val(strategy), val(model), path("assembly.fasta"), emit: assembly_fasta
    tuple val(sra), val(srr), path("assembly.gfa"), emit: assembly_graph
    tuple val(sra), val(srr), path("myloasm_*.log"), emit: assembly_log
    tuple val(sra), val(srr), path("assembly.bam"), path("assembly.bam.csi"), emit: assembly_bam

    script:
    """
    # Run myloasm
    myloasm ${reads} -o . -t ${task.cpus} --hifi

    # Rename outputs
    mv -v assembly_primary.fa assembly.fasta
    mv -v final_contig_graph.gfa assembly.gfa

    # Run minimap2
    minimap2 -ax map-hifi -t ${task.cpus} assembly.fasta ${reads} \\
      | samtools sort --output-fmt BAM -@ ${task.cpus} -o assembly.bam

    samtools index -c -o assembly.bam.csi -@ ${task.cpus} assembly.bam
    """
}


process DIAMOND {
    tag { "${sra}:${srr}" }
    publishDir "${params.outdir}/${sra}/${srr}/", mode: 'copy', overwrite: true

    input:
    tuple val(sra), val(srr), val(assembler), val(strategy), val(model), path(assembly_fasta)
    path uniprot_db

    output:
    tuple val(sra), val(srr), path("assembly_vs_uniprot.tsv"), emit: blast

    script:
    """
		diamond blastx --sensitive --query "${assembly_fasta}" \\
			--out "assembly_vs_uniprot.tsv" --db "${uniprot_db}" \\
			--outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \\
			--verbose --threads ${task.cpus} --evalue 1e-25 --max-target-seqs 5
    """
}


process BLOBTOOLS {
    tag { "${sra}:${srr}" }
    publishDir "${params.outdir}/${sra}/${srr}/", mode: 'copy', overwrite: true

    input:
    tuple val(sra), val(srr), val(assembler), val(strategy), val(model), path(assembly_fasta), path(blast), path(assembly_bam), path(assembly_csi)
    path taxdump

    output:
    tuple val(sra), val(srr), val(assembler), val(strategy), val(model), path(assembly_fasta), path("blobtools.csv"), emit: blobtable
    tuple val(sra), val(srr), path("blobtools*.svg"), optional:true, emit: blobplots

    script:
    """
		blobtools create --fasta "${assembly_fasta}" \
			--cov "${assembly_bam}" --hits "${blast}" --taxdump "${taxdump}" \
			--threads  ${task.cpus} 'blobtools'

    blobtools filter --table 'blobtools.tsv' \
      --table-fields gc,length,ncount,assembly_cov,bestsumorder_superkingdom,bestsumorder_kingdom,bestsumorder_phylum,bestsumorder_class,bestsumorder_order,bestsumorder_family,bestsumorder_genus,bestsumorder_species \
			'blobtools'

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
    tuple val(sra), val(srr), val(assembler), val(strategy), val(model), path(assembly_fasta), path(blobtable)
    path(taxa)
    path(taxdump)

    output:
    tuple val(sra), val(srr), val(assembler), val(strategy), val(model), path("summary.csv"), emit: summary
    tuple val(sra), val(srr), path("*.ids.csv"), optional:true,emit: extracted_ids
    tuple val(sra), val(srr), path("*.fasta"), optional:true, emit: extracted_fasta

    script:
    """
    extract_records.py --blobtable "${blobtable}" --fasta "${assembly_fasta}" \\
      --taxa "${taxa}" --taxdump "${taxdump}"
    """
}

    {
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
  sra_ch = Channel.fromPath(params.sra, checkIfExists: true)
                  .splitCsv(header: true, strip: true)
                  .map { row -> row.sra.trim() }
                  .filter { it }
                  .distinct()
  taxa_ch    = Channel.value( file(params.taxa) )
  uniprot_db_ch = Channel.value( file(params.uniprot_db) )
  taxdump_ch    = Channel.value( file(params.taxdump) )

  // TODO: add validate taxa
  validated_taxa = VALIDATE_TAXA(taxa_ch, taxdump_ch)

  // Step 1: extract metadata & filter SRR
  filtered_srr = DOWNLOAD_SRA_METADATA(sra_ch, validated_taxa)

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
  short_ch       = srr_reads.filter { sra, srr, platform, strategy, model, assembler, reads -> assembler.equalsIgnoreCase('short') }
  long_nano_ch   = srr_reads.filter { sra, srr, platform, strategy, model, assembler, reads -> assembler.equalsIgnoreCase('long_nano') }
  long_pacbio_ch = srr_reads.filter { sra, srr, platform, strategy, model, assembler, reads -> assembler.equalsIgnoreCase('long_pacbio') }
  long_hifi_ch   = srr_reads.filter { sra, srr, platform, strategy, model, assembler, reads -> assembler.equalsIgnoreCase('long_hifi') }

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
  fasta_by  = asm_fasta_ch.map   { sra, srr, asm, strategy, mdl, fasta -> tuple([sra,srr], [asm, strategy, mdl, fasta]) }
  blast_by  = diamond_ch.map     { sra, srr, hits  -> tuple([sra,srr], hits) }
  bam_by    = bam_ch.map         { sra, srr, bam, bai -> tuple([sra,srr], [bam, bai]) }

  // Join (fasta * diamond) then * bam
  fasta_blast = fasta_by.join(blast_by)
  fasta_blast_bam = fasta_blast.join(bam_by)
  // -> ( [sra,srr,asm], fasta, hits, [bam,bai] )

  // Unkey + call BlobTools
  blobtools_in = fasta_blast_bam.map { key, fasta, hits, pair ->
    def (sra, srr) = key
    def (asm, strategy, mdl, assembly) = fasta
    def (bam, bai) = pair
    tuple(sra, srr, asm, strategy, mdl, assembly, hits, bam, bai)
  }

  blobtools_result = BLOBTOOLS(blobtools_in, taxdump_ch)
  // blobtools_result.view { sra, srr, asm, tbl -> "BLOBTOOLS:\t${sra}\t${srr}\t${asm}\t${tbl}" }

  // Step 6: extract taxa
  extracted_taxa = EXTRACT_TAXA(blobtools_result.blobtable, validated_taxa, taxdump_ch)
}
