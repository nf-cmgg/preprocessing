#!/usr/bin/env nextflow

// This workflow is still missing proper handling of metrics files.
// Elprep filter supports generating intermediate metrics files (--mark-optical-duplicates-intermediate) which can be merged into a single file per sample.

nextflow.enable.dsl = 2

include { ELPREP_SPLIT  } from "../../modules/nf-core/modules/elprep/split"
include { ELPREP_FILTER } from "../../modules/nf-core/modules/elprep/filter"
include { ELPREP_MERGE  } from "../../modules/nf-core/modules/elprep/merge"

workflow bamprocessing_elprep {
    take:
        ch_bam_chunks

    main:
        ch_versions = Channel.empty()

        ELPREP_SPLIT(ch_bam_chunks).out.bam
        ch_versions = ch_versions.mix(ELPREP_SPLIT.out.versions)

        ELPREP_FILTER(ELPREP_SPLIT.out.bam)
        ch_versions = ch_versions.mix(ELPREP_FILTER.out.versions)

        ELPREP_MERGE(ELPREP_FILTER.out.bam)
        ch_versions = ch_versions.mix(ELPREP_MERGE.out.versions)

    emit:
        bam: ELPREP_MERGE.out.
        metrics: ELPREP_MERGE.out.metrics
        versions = ch_versions
}
