#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {BAMTOOLS_SPLIT as BAMPROCESSING_MERGE_SPLIT            } from "../../../modules/bamtools/split/main"
include {BIOBAMBAM_BAMMARKDUPLICATES2 as BAMPROCESSING_MARKDUP  } from "../../../modules/nf-core/modules/biobambam/bamsormadup/main"
include {SAMTOOLS_MERGE as BAMPROCESSING_MERGE                  } from "../../../modules/nf-core/modules/samtools/merge/main"

workflow MARKDUP_PARALLEL {
    take:
        ch_input_bam    // [meta, [bam1,bam2,...]]

    main:
        ch_versions = Channel.empty()

        // merge the bams and split per chromosome or interval
        BAMPROCESSING_MERGE_SPLIT(ch_input_bam)
        ch_versions = ch_versions.mix(BAMPROCESSING_MERGE_SPLIT.out.versions)

        // markduplicates
        BAMPROCESSING_MARKDUP(BAMPROCESSING_MERGE_SPLIT.out.bam)
        ch_versions = ch_versions.mix(BAMPROCESSING_MARKDUP.out.versions)

        // re-merge bam
        BAMPROCESSING_MERGE(BAMPROCESSING_MARKDUP.out.bam)
        ch_versions = ch_versions.mix(BAMPROCESSING_MERGE.out.versions)

        // TODO re-merge metrics
        ch_metrics_merged = Channel.empty()

    emit:
        bam: BAMPROCESSING_MERGE.out.bam    // [meta, bam]
        metrics: ch_metrics_merged          // [meta, metrics]
        versions = ch_versions              // versions
}
