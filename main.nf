include { SAMTOOLS_INDEX } from './modules/nf-core/samtools/index/main'
include { IVAR_TRIM } from './modules/nf-core/ivar/trim/main'

workflow {
    main:
    // Load bam alignments
    if (params.samplesheet) {
      bam_ch = channel.fromPath(params.samplesheet)
        | splitCsv(header: true)
        | map { row ->
            tuple([id: row.sample], [file(row.bam)])
          }
    } else if (params.bam) {
        bam_ch = channel.fromPath(params.bam, checkIfExists: true)
          | map { bamfile -> tuple([id: bamfile.baseName], bamfile) }
    } else {
        error "Please specify either --samplesheet samplesheet.csv or --bam 'data/*.bam'"
    }
    bam_ch | SAMTOOLS_INDEX
    bam_indexed_ch = bam_ch
    | join(SAMTOOLS_INDEX.out.index)

    primer_ch = channel.fromPath(params.primers, checkIfExists: true)

    ivar_trim_input = bam_indexed_ch | combine(primer_ch)

    IVAR_TRIM(
      bam_indexed_ch,
      ivar_trim_input | map {n -> n.get(3)}
    )

    IVAR_TRIM.out.bam | view
}
