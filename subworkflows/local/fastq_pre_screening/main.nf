include { SINGLEM } from '../../../modules/local/singlem'

/*
 * Pre-screen FASTQ samples with SingleM when target taxa screening is enabled.
 * The emitted reads channel keeps the assembly tuple shape in both branches.
 */

workflow FASTQ_PRE_SCREENING {
  take:
    fastq_samplesheet_channel
    validated_taxa_ch
    singlem_db_ch

  main:
    def doScreening = params.taxa != null

    if (doScreening) {
      fastq_for_singlem_channel = fastq_samplesheet_channel.map { sra, srr, platform, model, strategy, read_type, reads ->
        tuple(sra, srr, platform, model, strategy, read_type, 'RUN_SINGLEM', reads)
      }

      singlem_out = SINGLEM(fastq_for_singlem_channel, validated_taxa_ch, singlem_db_ch)
      screened_reads_channel = singlem_out.reads
      singlem_note_channel = singlem_out.note
    }
    else {
      screened_reads_channel = fastq_samplesheet_channel
      singlem_note_channel = channel.empty()
    }

  emit:
    reads = screened_reads_channel
    note = singlem_note_channel
}
