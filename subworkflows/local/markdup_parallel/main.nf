#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {BAMTOOLS_SPLIT               } from "../../../modules/nf-core/modules/bamtools/split/main"
include {BIOBAMBAM_BAMMARKDUPLICATES2 } from "../../../modules/nf-core/modules/biobambam/bammarkduplicates2/main"
include {SAMTOOLS_MERGE               } from "../../../modules/nf-core/modules/samtools/merge/main"

workflow MARKDUP_PARALLEL {
    take:
        ch_input_bam    // [meta, [bam1,bam2,...]]

    main:
        ch_versions = Channel.empty()

        // merge the bams and split per chromosome or interval
        BAMTOOLS_SPLIT(ch_input_bam)
        ch_versions = ch_versions.mix(BAMTOOLS_SPLIT.out.versions)

        // markduplicates
        BIOBAMBAM_BAMMARKDUPLICATES2(BAMTOOLS_SPLIT.out.bam)
        ch_versions = ch_versions.mix(BAMPROCESSING_MARKDUP.out.versions)

        // re-merge bam
        SAMTOOLS_MERGE(BIOBAMBAM_BAMMARKDUPLICATES2.out.bam, [])
        ch_versions = ch_versions.mix(BIOBAMBAM_BAMMARKDUPLICATES2.out.versions)

        // TODO re-merge metrics
        ch_metrics_merged = Channel.empty()

    emit:
        bam      = SAMTOOLS_MERGE.out.bam  // [meta, bam]
        metrics  = ch_metrics_merged       // [meta, metrics]
        versions = ch_versions             // versions
}
