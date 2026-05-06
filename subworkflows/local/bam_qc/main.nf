// samtools modules
include { SAMTOOLS_STATS                } from '../../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_IDXSTATS             } from '../../../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_FLAGSTAT             } from '../../../modules/nf-core/samtools/flagstat/main'

// picard modules
include { PICARD_COLLECTMULTIPLEMETRICS } from '../../../modules/nf-core/picard/collectmultiplemetrics/main'
include { PICARD_COLLECTHSMETRICS       } from '../../../modules/nf-core/picard/collecthsmetrics/main'
include { PICARD_COLLECTWGSMETRICS      } from '../../../modules/nf-core/picard/collectwgsmetrics/main'

workflow BAM_QC {
    take:
    ch_bam_bai_roi_fasta_fai_dict // channel: [ val(meta), path(bam), path(bai), path(roi), path(fasta), path(fai), path(dict)]

    main:
    ch_bam_bai_roi_fasta_fai_dict
        .map { meta, bam, bai, _roi, fasta, _fai, _dict ->
            return [meta, bam, bai, fasta]
        }
        .set { ch_bam_bai_fasta }

    SAMTOOLS_STATS(ch_bam_bai_fasta)
    SAMTOOLS_FLAGSTAT(ch_bam_bai_fasta.map { meta, bam, bai, _fasta -> return [meta, bam, bai] })
    SAMTOOLS_IDXSTATS(ch_bam_bai_fasta.map { meta, bam, bai, _fasta -> return [meta, bam, bai] })

    ch_picard_hsmetrics = channel.empty()
    ch_picard_multiplemetrics = channel.empty()
    ch_picard_multiplemetrics_pdf = channel.empty()
    ch_picard_wgsmetrics = channel.empty()

    ch_bam_bai_roi_fasta_fai_dict
        .filter { meta, _bam, _bai, _roi, _fasta, _fai, _dict ->
            !meta.disable_picard_metrics
        }
        .set { ch_picard }

    PICARD_COLLECTMULTIPLEMETRICS(ch_picard)
    ch_picard_multiplemetrics = PICARD_COLLECTMULTIPLEMETRICS.out.metrics
    ch_picard_multiplemetrics_pdf = PICARD_COLLECTMULTIPLEMETRICS.out.pdf

    ch_picard
        .branch { meta, bam, bai, roi, fasta, fai, dict ->
            hsmetrics: roi != []
            return [meta, bam, bai, roi, roi, fasta, fai, dict]
            wgsmetrics: roi == []
            return [meta, bam, bai, fasta, fai, dict]
        }
        .set { ch_picard_coverage }

    PICARD_COLLECTWGSMETRICS(ch_picard_coverage.wgsmetrics, [])
    ch_picard_wgsmetrics = PICARD_COLLECTWGSMETRICS.out.metrics

    PICARD_COLLECTHSMETRICS(ch_picard_coverage.hsmetrics)
    ch_picard_hsmetrics = PICARD_COLLECTHSMETRICS.out.metrics

    emit:
    samtools_stats             = SAMTOOLS_STATS.out.stats
    samtools_flagstat          = SAMTOOLS_FLAGSTAT.out.flagstat
    samtools_idxstats          = SAMTOOLS_IDXSTATS.out.idxstats
    picard_multiplemetrics     = ch_picard_multiplemetrics
    picard_multiplemetrics_pdf = ch_picard_multiplemetrics_pdf
    picard_wgsmetrics          = ch_picard_wgsmetrics
    picard_hsmetrics           = ch_picard_hsmetrics
}
