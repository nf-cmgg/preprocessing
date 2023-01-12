#!/usr/bin/env nextflow

include { MOSDEPTH } from "../../../modules/nf-core/mosdepth/main"

workflow COVERAGE {
    take:
        ch_reads_index      // channel: [mandatory][ meta,  reads, index ]
        ch_fasta_fai        // channel: [mandatory][ meta2, fasta, fai ]
        ch_target_interval  // channel: [optional] [ meta3, target_interval_bed ]

    main:
        ch_versions = Channel.empty()
        ch_metrics  = Channel.empty()

        ch_fasta    = ch_fasta_fai.map  {meta, fasta, fai -> [meta, fasta] }.collect()
        ch_meta_target_interval = ch_target_interval.map { it -> [[:],[]] }

        MOSDEPTH(ch_reads_index, ch_meta_target_interval, ch_fasta)
        ch_metrics = ch_metrics.mix(
            MOSDEPTH.out.summary_txt,
            MOSDEPTH.out.global_txt,
            MOSDEPTH.out.regions_txt
        ).groupTuple().map{ meta, metrics -> return [meta, metrics.flatten()]}

        ch_versions = ch_versions.mix(MOSDEPTH.out.versions)
        ch_metrics.dump(tag: "COVERAGE: metrics", {FormattingService.prettyFormat(it)})

    emit:
        per_base_bed        = MOSDEPTH.out.per_base_bed     // [meta, bed]
        per_base_bed_csi    = MOSDEPTH.out.per_base_csi     // [meta, csi]
        regions_bed         = MOSDEPTH.out.regions_bed      // [meta, bed]
        regions_bed_csi     = MOSDEPTH.out.regions_csi      // [meta, csi]
        quantized_bed       = MOSDEPTH.out.quantized_bed    // [meta, bed]
        quantized_bed_csi   = MOSDEPTH.out.quantized_csi    // [meta, csi]
        metrics             = ch_metrics                    // [[meta, [metrics]]
        versions            = ch_versions                   // [versions]
}
