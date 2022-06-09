#!/usr/bin/env nextflow

include { PICARD_COLLECTMULTIPLEMETRICS } from "../../../modules/nf-core/modules/picard/collectmultiplemetrics/main"
include { BAM_STATS_SAMTOOLS            } from "../../nf-core/subworkflows/bam_stats_samtools/main"

workflow BAM_QC {
    take:
        ch_bam_bai      // [meta, bam, bai]

    main:
        ch_versions = Channel.empty()
        ch_metrics  = Channel.empty()

        // Collect multiple metrics
        // PICARD_COLLECTMULTIPLEMETRICS( [meta, bam], fasta)
        PICARD_COLLECTMULTIPLEMETRICS( ch_bam_bai.map {meta, bam, bai -> return [meta,bam]}, [] )
        ch_metrics  = ch_metrics.mix(PICARD_COLLECTMULTIPLEMETRICS.out.metrics)
        ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions)

        // SUBWORKFLOW: bam_stats_samtools
        // Run samtools QC modules
        // BAM_STATS_SAMTOOLS([meta, bam, bai])
        BAM_STATS_SAMTOOLS(ch_bam_bai)
        ch_metrics = ch_metrics.mix(
            BAM_STATS_SAMTOOLS.out.stats,
            BAM_STATS_SAMTOOLS.out.flagstat,
            BAM_STATS_SAMTOOLS.out.idxstats
        )
        ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)

    emit:
        metrics  = ch_metrics  // [[meta, metrics], [...], ...]
        versions = ch_versions // [versions]
}
