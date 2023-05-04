#!/usr/bin/env nextflow

include { BAM_QC_PICARD             } from '../../nf-core/bam_qc_picard/main'
include { BAM_STATS_SAMTOOLS        } from "../../nf-core/bam_stats_samtools/main"

include { PICARD_BEDTOINTERVALLIST  } from '../../../modules/nf-core/picard/bedtointervallist/main'

workflow BAM_QC {
    take:
        ch_bam_bai_target   // channel: [mandatory] [meta, bam, bai, target]
        ch_fasta_fai        // channel: [mandatory] [meta2, fasta, fai]
        ch_fasta_dict       // channel: [mandatory] [meta2, dict]
        disable_picard      // boolean: [optional]  [true]

    main:
        ch_versions = Channel.empty()
        ch_metrics  = Channel.empty()

        ch_bait_interval_list = []
        ch_target_interval_list = []

        ch_fasta = ch_fasta_fai.map {meta, fasta, fai -> fasta             }.collect()
        ch_meta_fai   = ch_fasta_fai.map {meta, fasta, fai -> [meta, fai]  }.collect()
        ch_meta_fasta = ch_fasta_fai.map {meta, fasta, fai -> [meta, fasta]}.collect()

        if (!disable_picard) {

            PICARD_COLLECTMULTIPLEMETRICS( ch_bam_bai, ch_fasta, ch_fasta_fai )
            ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions.first())
            ch_metrics = ch_metrics.mix(PICARD_COLLECTMULTIPLEMETRICS.out.metrics)


            ch_bam_bai_target_branched = ch_bam_bai_target.branch {
                hsmetrics  : it.size == 4 && it[3] != []
                    return it
                wgsmetrics : true
                    return [ it[0], it[1], it[2] ]
            }

            // WGS metrics
            PICARD_COLLECTWGSMETRICS( ch_bam_bai_target_branched.wgsmetrics, ch_fasta, ch_fasta_fai, [] )
            ch_versions = ch_versions.mix(PICARD_COLLECTWGSMETRICS.out.versions.first())
            ch_metrics = ch_metrics.mix(PICARD_COLLECTWGSMETRICS.out.metrics)

            // HS metrics
            PICARD_COLLECTHSMETRICS( ch_bam_bai_target_branched.hsmetrics, ch_fasta, ch_fasta_fai, ch_fasta_dict)
            ch_coverage_metrics = ch_coverage_metrics.mix(PICARD_COLLECTHSMETRICS.out.metrics)
            ch_versions = ch_versions.mix(PICARD_COLLECTHSMETRICS.out.versions.first())


        }

        // SUBWORKFLOW: bam_stats_samtools
        // Run samtools QC modules
        // BAM_STATS_SAMTOOLS([meta, bam, bai])
        BAM_STATS_SAMTOOLS(ch_bam_bai, ch_fasta)

        ch_metrics = ch_metrics.mix(
            BAM_STATS_SAMTOOLS.out.stats,
            BAM_STATS_SAMTOOLS.out.flagstat,
            BAM_STATS_SAMTOOLS.out.idxstats
        )
        ch_metrics.dump(tag: "BAM QC: metrics", {FormattingService.prettyFormat(it)})
        ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)

    emit:
        metrics  = ch_metrics   // [[meta, metrics], [...], ...]
        versions = ch_versions  // [versions]
}
