#!/usr/bin/env nextflow

include { MOSDEPTH                  } from "../../../modules/nf-core/mosdepth/main"
include { MULTIQC as COVERAGE_MQC   } from "../../../modules/nf-core/multiqc/main"

include { PANELCOVERAGE             } from "../../../modules/local/panelcoverage/main"

workflow COVERAGE {
    take:
        ch_reads_index      // channel: [mandatory][ meta, reads, index ]
        ch_fasta_fai        // channel: [mandatory][ meta, fasta, fai ]
        ch_roi              // channel: [optional] [ meta, roi_bed ]
        ch_panels           // channel: [optional] [ panel_bed, ... ]

    main:
        ch_versions             = Channel.empty()
        ch_metrics              = Channel.empty()
        ch_coverage_mqc_files   = Channel.empty()

        ch_fasta      = ch_fasta_fai.map  {meta, fasta, fai -> [meta, fasta] }.collect()

        // Get per base coverage
        MOSDEPTH(ch_reads_index, ch_roi, ch_fasta)
        ch_metrics = ch_metrics
            .mix(
                MOSDEPTH.out.summary_txt,
                MOSDEPTH.out.global_txt,
                MOSDEPTH.out.regions_txt
            )
            .groupTuple()
            .map { meta, metrics -> [meta, metrics.flatten()] }

        ch_versions = ch_versions.mix(MOSDEPTH.out.versions)
        ch_metrics.dump(tag: "COVERAGE: metrics", {FormattingService.prettyFormat(it)})

        if (ch_panels) {
            // Mock Mosdepth distribution output
            PANELCOVERAGE(MOSDEPTH.out.per_base_bed, ch_panels)
            ch_versions = ch_versions.mix(PANELCOVERAGE.out.versions)
            ch_coverage_mqc_files = ch_coverage_mqc_files.mix(PANELCOVERAGE.out.mqc_files).collect()

            // Generate Coverage MultiQC
            COVERAGE_MQC(ch_coverage_mqc_files)
            ch_versions = ch_versions.mix(MULTIQC.out.versions)
        }


    emit:
        per_base_bed        = MOSDEPTH.out.per_base_bed         // [meta, bed]
        per_base_bed_csi    = MOSDEPTH.out.per_base_csi         // [meta, csi]
        regions_bed         = MOSDEPTH.out.regions_bed          // [meta, bed]
        regions_bed_csi     = MOSDEPTH.out.regions_csi          // [meta, csi]
        quantized_bed       = MOSDEPTH.out.quantized_bed        // [meta, bed]
        quantized_bed_csi   = MOSDEPTH.out.quantized_csi        // [meta, csi]
        metrics             = ch_metrics                        // [[meta, [metrics]]
        versions            = ch_versions                       // [versions]
}
