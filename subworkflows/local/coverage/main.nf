#!/usr/bin/env nextflow

include { MOSDEPTH                   } from "../../../modules/nf-core/modules/mosdepth/main"
include { PICARD_COLLECTWGSMETRICS   } from "../../../modules/nf-core/modules/picard/collectwgsmetrics/main"
include { PICARD_COLLECTHSMETRICS    } from "../../../modules/nf-core/modules/picard/collecthsmetrics/main"

workflow COVERAGE {
    take:
        ch_bam_bai
        ch_fasta
        ch_fai
        ch_target_interval
        ch_bait_interval

    main:
        ch_versions = Channel.empty()
        ch_metrics  = Channel.empty()

        MOSDEPTH(ch_bam_bai, ch_target_interval, ch_fasta)
        ch_metrics = ch_metrics.mix(
            MOSDEPTH.out.summary_txt,
            MOSDEPTH.out.global_txt,
            MOSDEPTH.out.regions_txt
        )

        ch_bam = ch_bam_bai.map { meta, bam, bai -> return [meta, bam]}
        if (ch_bait_interval || ch_target_interval) {
            if (!ch_bait_interval) log.error("Bait interval channel is empty")
            if (!ch_target_interval) log.error("Target interval channel is empty")
            PICARD_COLLECTHSMETRICS( ch_bam, ch_fasta, ch_fasta_fai, ch_bait_interval, ch_target_interval )
            ch_metrics = ch_metrics.mix(PICARD_COLLECTHSMETRICS.out.metrics)
            ch_versions = ch_versions.mix(PICARD_COLLECTHSMETRICS.out.versions)
        } else {
            PICARD_COLLECTWGSMETRICS( ch_bam, ch_fasta )
            ch_versions = ch_versions.mix(PICARD_COLLECTWGSMETRICS.out.versions)
            ch_metrics = ch_metrics.mix(PICARD_COLLECTWGSMETRICS.out.metrics)
        }

    emit:
        per_base_bed        = MOSDEPTH.out.per_base_bed  // [meta, bed]
        per_base_bed_csi    = MOSDEPTH.out.per_base_csi  // [meta, csi]
        regions_bed         = MOSDEPTH.out.regions_bed   // [meta, bed]
        regions_bed_csi     = MOSDEPTH.out.regions_csi   // [meta, csi]
        quantized_bed       = MOSDEPTH.out.quantized_bed // [meta, bed]
        quantized_bed_csi   = MOSDEPTH.out.quantized_csi // [meta, csi]
        metrics             = ch_metrics                 // [[meta, metrics],...]
        versions            = ch_versions                // [versions]
}
