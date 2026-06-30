include { MOSDEPTH          } from "../../../modules/nf-core/mosdepth/main.nf"
include { PANELCOVERAGE     } from "../../../modules/local/panelcoverage/main"
include { RIKER_MULTI       } from '../../../modules/nf-core/riker/multi/main'
include { SAMTOOLS_COVERAGE } from "../../../modules/nf-core/samtools/coverage/main"
include { SAMTOOLS_FLAGSTAT } from '../../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_IDXSTATS } from '../../../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_STATS    } from '../../../modules/nf-core/samtools/stats/main'

workflow BAM_QC {
    take:
    ch_bam_bai_roi_fasta_fai // channel: [ val(meta), path(bam), path(bai), path(roi), path(fasta), path(fai)]
    ch_genelists // channel: [optional] [genelists]

    main:

    ch_bam_bai_roi_fasta_fai
        .map { meta, bam, bai, _roi, fasta, fai ->
            return [meta, bam, bai, fasta, fai]
        }
        .set { ch_bam_bai_fasta_fai }

    // basic QC
    SAMTOOLS_FLAGSTAT(ch_bam_bai_fasta_fai.map { meta, bam, bai, _fasta, _fai ->
        return [meta, bam, bai]
    })
    SAMTOOLS_IDXSTATS(ch_bam_bai_fasta_fai.map { meta, bam, bai, _fasta, _fai ->
        return [meta, bam, bai]
    })

    MOSDEPTH(
        ch_bam_bai_roi_fasta_fai.map { meta, bam, bai, roi, fasta, _fai ->
            return [meta, bam, bai, roi, fasta]
        },
        ['NO_COVERAGE', 'LOW_COVERAGE', 'CALLABLE']
    )

    // full QC
    // Only run on samples requiring full QC
    ch_full_qc = ch_bam_bai_roi_fasta_fai.filter { meta, _bam, _bai, _roi, _fasta, _fai ->
        meta.qc_mode == "full"
    }

    SAMTOOLS_STATS(ch_full_qc.map { meta, bam, bai, _roi, fasta, fai ->
        return [meta, bam, bai, fasta, fai]
    })

    RIKER_MULTI(ch_full_qc)

    SAMTOOLS_COVERAGE(
        ch_full_qc.map { meta, bam, bai, _roi, fasta, fai ->
            return [meta, bam, bai, fasta, fai]
        }
    )

    PANELCOVERAGE(
        MOSDEPTH.out.per_base_bed.join(MOSDEPTH.out.per_base_csi)
        .combine(ch_genelists)
        .filter{ meta, _bed, _index, _genelists -> meta.qc_mode == "full" }
        .map { meta, bed, index, genelists ->
            // Because groovy typing sucks ass; apparently an array of 1 is automatically converted to a string...
            if (genelists !instanceof List) {
                genelists = [genelists]
            }
            def filtered_genelists = (meta.tag && meta.tag.toLowerCase() == "seqcap")
                ? genelists.findAll { genelist -> genelist.name.toLowerCase().contains("seqcap") }
                : genelists.findAll { genelist -> !genelist.name.toLowerCase().contains("seqcap") }

            if (filtered_genelists.size() > 0) {
                return [
                    meta,
                    bed,
                    index,
                    filtered_genelists,
                ]
            }
        }
    )


    emit:
    mosdepth_global         = MOSDEPTH.out.global_txt
    mosdepth_per_base_bed   = MOSDEPTH.out.per_base_bed
    mosdepth_per_base_csi   = MOSDEPTH.out.per_base_csi
    mosdepth_per_base_d4    = MOSDEPTH.out.per_base_d4
    mosdepth_quantized_bed  = MOSDEPTH.out.quantized_bed
    mosdepth_quantized_csi  = MOSDEPTH.out.quantized_csi
    mosdepth_regions        = MOSDEPTH.out.regions_txt
    mosdepth_regions_bed    = MOSDEPTH.out.regions_bed
    mosdepth_regions_csi    = MOSDEPTH.out.regions_csi
    mosdepth_summary        = MOSDEPTH.out.summary_txt
    mosdepth_thresholds_bed = MOSDEPTH.out.thresholds_bed
    mosdepth_thresholds_csi = MOSDEPTH.out.thresholds_csi
    panelcoverage           = PANELCOVERAGE.out.regiondist
    riker_metrics           = RIKER_MULTI.out.metrics
    samtools_coverage       = SAMTOOLS_COVERAGE.out.coverage
    samtools_flagstat       = SAMTOOLS_FLAGSTAT.out.flagstat
    samtools_idxstats       = SAMTOOLS_IDXSTATS.out.idxstats
    samtools_stats          = SAMTOOLS_STATS.out.stats
}
