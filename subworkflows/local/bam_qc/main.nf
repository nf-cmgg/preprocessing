// samtools modules
include { SAMTOOLS_STATS    } from '../../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_IDXSTATS } from '../../../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_FLAGSTAT } from '../../../modules/nf-core/samtools/flagstat/main'

// riker modules
include { RIKER_MULTI       } from '../../../modules/nf-core/riker/multi/main'

workflow BAM_QC {
    take:
    ch_bam_bai_roi_fasta_fai // channel: [ val(meta), path(bam), path(bai), path(roi), path(fasta), path(fai)]

    main:
    ch_bam_bai_roi_fasta_fai
        .map { meta, bam, bai, _roi, fasta, fai ->
            return [meta, bam, bai, fasta, fai]
        }
        .set { ch_bam_bai_fasta_fai }

    SAMTOOLS_STATS(ch_bam_bai_fasta_fai)
    SAMTOOLS_FLAGSTAT(ch_bam_bai_fasta_fai.map { meta, bam, bai, _fasta, _fai ->
        return [meta, bam, bai]
    })
    SAMTOOLS_IDXSTATS(ch_bam_bai_fasta_fai.map { meta, bam, bai, _fasta, _fai ->
        return [meta, bam, bai]
    })

    RIKER_MULTI(ch_bam_bai_roi_fasta_fai)

    emit:
    samtools_stats          = SAMTOOLS_STATS.out.stats
    samtools_flagstat       = SAMTOOLS_FLAGSTAT.out.flagstat
    samtools_idxstats       = SAMTOOLS_IDXSTATS.out.idxstats
    riker_metrics           = RIKER_MULTI.out.metrics
}
