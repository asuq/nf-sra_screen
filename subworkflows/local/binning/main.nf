include {
  validateBinningOptions
  binnersForSample
  binnerCsvContains
  binnerCsvSize
} from '../../../lib/workflow_helpers.nf'

include { METABAT } from '../../../modules/local/metabat'
include { COMEBIN } from '../../../modules/local/comebin'
include { COMEBIN_GPU } from '../../../modules/local/comebin_gpu'
include { VAMB } from '../../../modules/local/vamb'
include { VAMB_GPU } from '../../../modules/local/vamb_gpu'
include { LORBIN } from '../../../modules/local/lorbin'
include { LORBIN_GPU } from '../../../modules/local/lorbin_gpu'
include { SEMIBIN } from '../../../modules/local/semibin'
include { ROSELLA } from '../../../modules/local/rosella'
include { DASTOOL } from '../../../modules/local/dastool'
include { BINETTE } from '../../../modules/local/binette'

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
      .groupTuple(remainder: true)
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
