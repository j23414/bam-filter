include { SAMTOOLS_INDEX} from './modules/nf-core/samtools/index/main'
include { IVAR_TRIM } from './modules/nf-core/ivar/trim/main'

process SAMTOOLS_FILTER {
    tag "${meta.id}"
    label 'process_low'

    conda "./modules/nf-core/samtools/view/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

  input: tuple val(meta), path(bam)
  output:
  tuple val(meta), path("${bam.baseName}_filtered.bam"), emit: bam

  script:
  """
  samtools sort -n --threads ${task.cpus - 1} ${bam} \
    | samtools fixmate -m -@ 2 - - \
    | samtools sort --write-index -@ 2 - \
    | samtools markdup --write-index -@ 2 - ${bam.baseName}_filtered.bam
  """
}

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

    IVAR_TRIM.out.bam | SAMTOOLS_FILTER | view
}
