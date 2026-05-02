include { DOWNLOAD_SRA_METADATA } from '../../../modules/local/download_sra_metadata'
include { SANDPIPER } from '../../../modules/local/sandpiper'
include { DOWNLOAD_SRR } from '../../../modules/local/download_srr'
include { SINGLEM } from '../../../modules/local/singlem'

workflow PRE_SCREENING {
  take:
    sra_accessions_channel
    fastq_samplesheet_channel
    validated_taxa_ch
    sandpiper_db_ch
    singlem_db_ch

  main:
    def doScreening = params.taxa != null

    // Start PRE_SCREENING after taxa validation
    def gated_sra_accessions_channel = doScreening
                        ? sra_accessions_channel.combine(validated_taxa_ch).map { sra, vt -> sra }
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
      sandpiper = SANDPIPER(srr_ch, validated_taxa_ch, sandpiper_db_ch)
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

      sra_reads_for_singlem_channel = download_srr.reads.map { sra, srr, platform, model, strategy, detected_read_type, sandpiper_dec, reads, asm_txt ->
        def fixedReadType = file(asm_txt)?.text?.trim() ?: detected_read_type
        tuple(sra, srr, platform, model, strategy, fixedReadType, sandpiper_dec, reads)
      }

      fastq_reads_for_singlem_channel = fastq_samplesheet_channel.map { sra, srr, platform, model, strategy, read_type, reads ->
        tuple(sra, srr, platform, model, strategy, read_type, 'RUN_SINGLEM', reads)
      }

      singlem_input_channel = channel.empty()
        .mix(sra_reads_for_singlem_channel)
        .mix(fastq_reads_for_singlem_channel)

      // SINGLEM prescreening for all downloaded SRA reads and user-provided FASTQ rows.
      singlem = SINGLEM(singlem_input_channel, validated_taxa_ch, singlem_db_ch)
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
      sra_downloaded_reads_channel = download_srr.reads.map { sra, srr, platform, model, strategy, detected_read_type, sandpiper_dec, reads, asm_txt ->
        def fixedReadType = file(asm_txt)?.text?.trim() ?: detected_read_type
        tuple(sra, srr, platform, model, strategy, fixedReadType, reads)
      }
      singlem_reads_ch = channel.empty()
        .mix(sra_downloaded_reads_channel)
        .mix(fastq_samplesheet_channel)

      sandpiper_note_ch    = channel.empty()
      download_srr_note_ch = download_srr.note
      singlem_note_ch      = channel.empty()
    }

  emit:
    // For assembly
    reads                 = singlem_reads_ch

    // For summary
    sra_metadata_skipped  = sra_metadata.skipped_sra
    sra_metadata_note     = sra_metadata.note
    sandpiper_note        = sandpiper_note_ch
    download_srr_note     = download_srr_note_ch
    singlem_note          = singlem_note_ch
}
