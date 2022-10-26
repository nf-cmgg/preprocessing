#!/usr/bin/env nextflow

include { MOSDEPTH                   } from "../../../modules/nf-core/mosdepth/main"
include { PICARD_COLLECTWGSMETRICS   } from "../../../modules/nf-core/picard/collectwgsmetrics/main"
include { PICARD_COLLECTHSMETRICS    } from "../../../modules/nf-core/picard/collecthsmetrics/main"

workflow COVERAGE {
    take:
        ch_reads_index      // channel: [mandatory][ meta, reads, index ]
        ch_fasta_fai        // channel: [mandatory][ meta2, fasta, fai ]
        ch_fasta_dict       // channel: [mandatory][ meta2, dict ]
        ch_target_interval  // channel: [optional] [ target_interval_bed ]
        ch_bait_interval    // channel: [optional] [ bai_interval_bed ]

    main:
        ch_versions = Channel.empty()
        ch_metrics  = Channel.empty()

        ch_fai        = ch_fasta_fai.map {meta, fasta, fai -> fai          }.collect()
        ch_fasta      = ch_fasta_fai.map {meta, fasta, fai -> fasta        }.collect()
        ch_dict       = ch_fasta_dict.map {meta, dict      -> dict         }.collect()

        MOSDEPTH(ch_reads_index, ch_target_interval, ch_fasta)
        ch_metrics = ch_metrics.mix(
            MOSDEPTH.out.summary_txt,
            MOSDEPTH.out.global_txt,
            MOSDEPTH.out.regions_txt
        )
        ch_versions = ch_versions.mix(MOSDEPTH.out.versions)

        ch_reads = ch_reads_index.map { meta, reads, index -> return [meta, reads]}
        if (ch_bait_interval || ch_target_interval) {
            if (!ch_bait_interval) log.error("Bait interval channel is empty")
            if (!ch_target_interval) log.error("Target interval channel is empty")
            PICARD_COLLECTHSMETRICS( ch_reads, ch_fasta, ch_fai, ch_dict, ch_bait_interval, ch_target_interval )
            ch_metrics = ch_metrics.mix(PICARD_COLLECTHSMETRICS.out.metrics)
            ch_versions = ch_versions.mix(PICARD_COLLECTHSMETRICS.out.versions)
        } else {
            PICARD_COLLECTWGSMETRICS( ch_reads, ch_fasta )
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
