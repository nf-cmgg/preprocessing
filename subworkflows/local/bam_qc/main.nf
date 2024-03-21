// samtools modules
include { SAMTOOLS_STATS    } from '../../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_IDXSTATS } from '../../../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_FLAGSTAT } from '../../../modules/nf-core/samtools/flagstat/main'

// picard modules
include { PICARD_COLLECTMULTIPLEMETRICS } from '../../../modules/nf-core/picard/collectmultiplemetrics/main'
include { PICARD_COLLECTHSMETRICS       } from '../../../modules/nf-core/picard/collecthsmetrics/main'
include { PICARD_COLLECTWGSMETRICS      } from '../../../modules/nf-core/picard/collectwgsmetrics/main'

workflow BAM_QC {
    take:
    ch_bam_bai_roi_fasta_fai_dict   // channel: [ val(meta), path(bam), path(bai), path(roi), path(fasta), path(fai), path(dict)]
    disable_picard                  // boolean

    main:
    ch_versions = Channel.empty()

    ch_bam_bai_roi_fasta_fai_dict
    .map{ meta, bam, bai, roi, fasta, fai, dict -> return [meta, bam, bai, fasta]}
    .set{ ch_bam_bai_fasta }

    SAMTOOLS_STATS ( ch_bam_bai_fasta )
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions)

    ch_bam_bai_fasta
    .map{ meta, bam, bai, fasta -> return [meta, bam, bai]}
    .set{ ch_bam_bai }

    SAMTOOLS_FLAGSTAT ( ch_bam_bai )
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions)

    SAMTOOLS_IDXSTATS ( ch_bam_bai )
    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions)

    ch_picard_hsmetrics = Channel.empty()
    ch_picard_multiplemetrics = Channel.empty()
    ch_picard_wgsmetrics = Channel.empty()
    if (!disable_picard) {

        ch_bam_bai_roi_fasta_fai_dict
        .map{ meta, bam, bai, roi, fasta, fai, dict -> return [meta, bam, bai, fasta, fai]}
        .set{ ch_bam_bai_fasta_fai }

        PICARD_COLLECTMULTIPLEMETRICS ( ch_bam_bai_fasta_fai )
        ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions)
        ch_picard_multiplemetrics = ch_picard_multiplemetrics.mix(PICARD_COLLECTMULTIPLEMETRICS.out.metrics)

        ch_bam_bai_roi_fasta_fai_dict
        .branch{ meta, bam, bai, roi, fasta, fai, dict ->
            hsmetrics   : roi != []
                return [meta, bam, bai, roi, roi, fasta, fai, dict]
            wgsmetrics  : roi == []
                return [meta, bam, bai, fasta, fai]
        }
        .set{ch_picard}

        PICARD_COLLECTWGSMETRICS ( ch_picard.wgsmetrics, [] )
        ch_versions = ch_versions.mix(PICARD_COLLECTWGSMETRICS.out.versions)
        ch_picard_wgsmetrics = ch_picard_wgsmetrics.mix(PICARD_COLLECTWGSMETRICS.out.metrics)

        PICARD_COLLECTHSMETRICS ( ch_picard.hsmetrics )
        ch_versions = ch_versions.mix(PICARD_COLLECTHSMETRICS.out.versions)
        ch_picard_hsmetrics = ch_picard_hsmetrics.mix(PICARD_COLLECTHSMETRICS.out.metrics)
    }

    emit:
    samtools_stats          = SAMTOOLS_STATS.out.stats
    samtools_flagstat       = SAMTOOLS_FLAGSTAT.out.flagstat
    samtools_idxstats       = SAMTOOLS_IDXSTATS.out.idxstats
    picard_multiplemetrics  = ch_picard_multiplemetrics
    picard_wgsmetrics       = ch_picard_wgsmetrics
    picard_hsmetrics        = ch_picard_hsmetrics
    versions                = ch_versions
}
