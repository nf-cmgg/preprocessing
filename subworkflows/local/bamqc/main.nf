//
// Run QC steps on BAM files
//

params.options = [:]

include { SAMTOOLS_QC                   } from '../../modules/local/samtoolsqc/main'                                addParams( options: params.options )
include { PICARD_COLLECTMULTIPLEMETRICS } from '../../modules/nf-core/modules/picardcollectmultiplemetrics/main'    addParams( options: params.options )
include { PICARD_COLLECTWGSMETRICS      } from '../../modules/nf-core/modules/picardcollectwgsmetrics/main'         addParams( options: params.options )
include { PICARD_COLLECTHSMETRICS       } from '../../modules/nf-core/modules/picardcollecthsmetrics/main'          addParams( options: params.options )

workflow BAM_QC {
    take:
    ch_bam_bai          // channel: [ val(meta), [ bam ], [bai/csi] ]
    ch_fasta            // channel: [ fasta ]
    ch_bait_interval    // channel: [ bait_interval ]
    ch_target_interval  // channel: [ target_interval ]

    main:
    ch_versions = Channel.empty()

    SAMTOOLS_QC( ch_bam_bai, ch_fasta )
    ch_versions = ch_versions.mix(SAMTOOLS_QC.out.versions.first())

    PICARD_COLLECTMULTIPLEMETRICS( ch_bam_bai, ch_fasta] )
    ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions.first())
    if (params.targetted) {
        if (ch_bait_interval.isEmpty()) {
            throw new Error("Bait interval channel is empty")
        }
        if (ch_target_interval.isEmpty()) {
            throw new Error("Target interval channel is empty")
        }
        PICARD_COLLECTHSMETRICS( ch_bam_bai, ch_fasta, ch_bait_interval, ch_target_interval )
        ch_versions = ch_versions.mix(PICARD_COLLECTHSMETRICS.out.versions.first())
    } else {
        PICARD_COLLECTWGSMETRICS( ch_bam_bai, ch_fasta )
        ch_versions = ch_versions.mix(PICARD_COLLECTWGSMETRICS.out.versions.first())
    }

    emit:
    stats               = SAMTOOLS_QC.out.stats                     // channel: [ val(meta), [ stats ] ]
    flagstat            = SAMTOOLS_QC.out.flagstat                  // channel: [ val(meta), [ flagstat ] ]
    idxstats            = SAMTOOLS_QC.out.idxstats                  // channel: [ val(meta), [ idxstats ] ]
    hs_metrics          = PICARD_COLLECTHSMETRICS.out.hs_metrics    // channel: [ val(meta), [ hs_metrics ] ]
    wgs_metrics         = PICARD_COLLECTWGSMETRICS.out.wgs_metrics  // channel: [ val(meta), [ wgs_metrics ] ]
    multiple_metrics    = PICARD_COLLECTMULTIPLEMETRICS.out.metrics // channel: [ val(meta), [ multiple_metrics ] ]

    versions = ch_versions                                          // channel: [ versions.yml ]
}
