#!/usr/bin/env nextflow

include { MOSDEPTH              } from "../../../modules/nf-core/mosdepth/main"
include { PANEL_COVERAGE        } from "../../../modules/local/panel_coverage/main"
include { MULTIQC as COVERAGEQC } from "../../../modules/nf-core/multiqc/main"

workflow COVERAGE {
    take:
        ch_reads_index_target   // channel: [mandatory][ meta,  reads, index, target_bed ]
        ch_fasta_fai            // channel: [mandatory][ meta2, fasta, fai ]
        ch_genelists            // channel: [optional] [ genelist ]

    main:
        ch_versions = Channel.empty()
        ch_metrics  = Channel.empty()

        ch_meta_fasta    = ch_fasta_fai.map  {meta, fasta, fai -> [meta, fasta] }.collect()

        MOSDEPTH(ch_reads_index_target, ch_meta_fasta)
        ch_metrics = ch_metrics.mix(
            MOSDEPTH.out.summary_txt,
            MOSDEPTH.out.global_txt,
            MOSDEPTH.out.regions_txt
        )

        ch_versions = ch_versions.mix(MOSDEPTH.out.versions)
        ch_metrics.dump(tag: "COVERAGE: metrics", {FormattingService.prettyFormat(it)})

        // separate WES/WGS samples to run genelist coverage on them
        ch_per_base_bed = MOSDEPTH.out.per_base_bed.branch{
            genelist_coverage: ["WES", "WGS"].contains(it[0].tag)
            other: true
        }

        ch_per_base_genelist = ch_per_base_bed.genelist_coverage.combine(ch_genelists)
        ch_per_base_genelist.dump(tag: "COVERAGE: per base bed with genelist", {FormattingService.prettyFormat(it)})

        //PANEL_COVERAGE(per_base_genelist)
        PANEL_COVERAGE(ch_per_base_genelist)

        ch_regions_dist = PANEL_COVERAGE.out.region_dist.map{ meta, files -> files}.collect()

        COVERAGEQC(ch_regions_dist, [], [], [])

    emit:
        per_base_bed        = MOSDEPTH.out.per_base_bed     // [meta, bed]
        per_base_bed_csi    = MOSDEPTH.out.per_base_csi     // [meta, csi]
        regions_bed         = MOSDEPTH.out.regions_bed      // [meta, bed]
        regions_bed_csi     = MOSDEPTH.out.regions_csi      // [meta, csi]
        regions_dist        = PANEL_COVERAGE.out.region_dist// [meta, [dist,dist,...]]
        quantized_bed       = MOSDEPTH.out.quantized_bed    // [meta, bed]
        quantized_bed_csi   = MOSDEPTH.out.quantized_csi    // [meta, csi]
        metrics             = ch_metrics                    // [[meta, [metrics]]
        versions            = ch_versions                   // [versions]
}

