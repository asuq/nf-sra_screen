include { CREATE_EMPTY_SUMMARY } from '../../../modules/local/create_empty_summary'
include { APPEND_SUMMARY } from '../../../modules/local/append_summary'

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


    // 1) Convert skipped SRA metadata into per-SRR rows with a textual note
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
        // Empty note or GTDB-only "soft skip" => success
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

        def note_texts = note_entries
          .sort { left, right -> left[0] <=> right[0] }
          .collect { note_entry ->
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
      successful_summary_by_binning_key = summary_rows_with_taxa_notes_channel.map { sra, srr, platform, model, strategy, read_type, assembler, summary_csv, base_note ->
        def sample_key_map = [sra: sra, srr: srr, platform: platform,
                              model: model, strategy: strategy, read_type: read_type, assembler: assembler]
        tuple(sample_key_map, [summary_csv, base_note])
      }

      binning_note_by_success_key = binning_notes_by_sample_channel.map { sra, srr, platform, model, strategy, read_type, assembler, binning_note ->
        def sample_key_map = [sra: sra, srr: srr, platform: platform,
                              model: model, strategy: strategy, read_type: read_type, assembler: assembler]
        tuple(sample_key_map, binning_note)
      }

      success_and_binning_notes_channel = channel.empty()
        .mix(successful_summary_by_binning_key)
        .mix(binning_note_by_success_key)
        .groupTuple()

      final_success_rows_channel = success_and_binning_notes_channel
        .map { grouped_sample_key, grouped_payloads ->
          def sample_key = grouped_sample_key as Map
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
