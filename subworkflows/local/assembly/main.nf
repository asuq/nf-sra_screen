include { assemblersForReadType } from '../../../lib/workflow_helpers.nf'

include { FASTP_SHORT } from '../../../modules/local/fastp_short'
include { METASPADES } from '../../../modules/local/metaspades'
include { UNICYCLER } from '../../../modules/local/unicycler'
include { METAFLYE_NANO } from '../../../modules/local/metaflye_nano'
include { METAFLYE_PACBIO } from '../../../modules/local/metaflye_pacbio'
include { METAFLYE_HIFI } from '../../../modules/local/metaflye_hifi'
include { MYLOASM } from '../../../modules/local/myloasm'
include { MAP_TO_ASSEMBLY } from '../../../modules/local/map_to_assembly'
include { DIAMOND } from '../../../modules/local/diamond'
include { BLOBTOOLS } from '../../../modules/local/blobtools'
include { EXTRACT_TAXA } from '../../../modules/local/extract_taxa'
include { CREATE_EMPTY_SUMMARY } from '../../../modules/local/create_empty_summary'
include { CREATE_ASSEMBLER_SELECTION_NOTE } from '../../../modules/local/create_assembler_selection_note'

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

