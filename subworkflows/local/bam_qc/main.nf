#!/usr/bin/env nextflow

include { BAM_QC_PICARD         } from '../../nf-core/bam_qc_picard/main'
include { BAM_STATS_SAMTOOLS    } from "../../nf-core/bam_stats_samtools/main"

include { PICARD_BEDTOINTERVALLIST as BAITTOINTERVALLIST    } from '../../../modules/nf-core/picard/bedtointervallist/main'
include { PICARD_BEDTOINTERVALLIST as TARGETTOINTERVALLIST  } from '../../../modules/nf-core/picard/bedtointervallist/main'

workflow BAM_QC {
    take:
        ch_bam_bai          // channel: [mandatory] [meta, bam, bai]
        ch_fasta_fai        // channel: [mandatory] [meta2, fasta, fai]
        ch_fasta_dict       // channel: [mandatory][ meta2, dict ]
        ch_target_interval  // channel: [optional] [ target_interval_bed ]
        ch_bait_interval    // channel: [optional] [ bai_interval_bed ]
        disable_picard      // boolean: [optional] [ true ]

    main:
        ch_versions = Channel.empty()
        ch_metrics  = Channel.empty()

        ch_bait_interval_list = []
        ch_target_interval_list = []

        ch_fasta = ch_fasta_fai.map {meta, fasta, fai -> fasta             }.collect()
        ch_meta_fai   = ch_fasta_fai.map {meta, fasta, fai -> [meta, fai]  }.collect()
        ch_meta_fasta = ch_fasta_fai.map {meta, fasta, fai -> [meta, fasta]}.collect()

        if (!disable_picard) {
            if ( ch_bait_interval ) {
                BAITTOINTERVALLIST(
                    ch_bait_interval.map{it -> return [[id:"bait"], it ]},
                    ch_fasta_dict,
                    []
                )
                ch_versions = ch_versions.mix(BAITTOINTERVALLIST.out.versions)
                ch_bait_interval_list = BAITTOINTERVALLIST.out.interval_list.map{meta, list -> list}.collect()
            }
            if ( ch_target_interval ) {
                TARGETTOINTERVALLIST(
                    ch_target_interval.map{[[id:"target"], it ]},
                    ch_fasta_dict,
                    []
                )
                ch_versions = ch_versions.mix(TARGETTOINTERVALLIST.out.versions)
                ch_target_interval_list = TARGETTOINTERVALLIST.out.interval_list.map{meta, list -> list}.collect()
            }

            // SUBWORKFLOW: bam_qc_picard
            // Run Picard QC modules
            // BAM_QC_PICARD([meta, bam, bai], [meta2, fasta], [meta2, fai], [bait_interval], [target_interval])
            BAM_QC_PICARD( ch_bam_bai, ch_meta_fasta, ch_meta_fai, ch_bait_interval_list, ch_target_interval_list )
            ch_metrics  = ch_metrics.mix(BAM_QC_PICARD.out.coverage_metrics, BAM_QC_PICARD.out.multiple_metrics)
            ch_versions = ch_versions.mix(BAM_QC_PICARD.out.versions)
        }

        // SUBWORKFLOW: bam_stats_samtools
        // Run samtools QC modules
        // BAM_STATS_SAMTOOLS([meta, bam, bai])
        BAM_STATS_SAMTOOLS(ch_bam_bai, ch_fasta)

        ch_metrics = ch_metrics.mix(
            BAM_STATS_SAMTOOLS.out.stats,
            BAM_STATS_SAMTOOLS.out.flagstat,
            BAM_STATS_SAMTOOLS.out.idxstats
        ).groupTuple().map { meta, metrics -> [meta, metrics.flatten()] }
        ch_metrics.dump(tag: "BAM QC: metrics", {FormattingService.prettyFormat(it)})
        ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)

    emit:
        metrics  = ch_metrics   // [[meta, metrics], [...], ...]
        versions = ch_versions  // [versions]
}
