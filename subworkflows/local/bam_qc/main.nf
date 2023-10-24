#!/usr/bin/env nextflow

include { BAM_STATS_SAMTOOLS        } from "../../nf-core/bam_stats_samtools/main"

include { PICARD_COLLECTMULTIPLEMETRICS } from '../../../modules/nf-core/picard/collectmultiplemetrics/main'
include { PICARD_COLLECTHSMETRICS       } from '../../../modules/nf-core/picard/collecthsmetrics/main'
include { PICARD_COLLECTWGSMETRICS      } from '../../../modules/nf-core/picard/collectwgsmetrics/main'


workflow BAM_QC {
    take:
        ch_meta_reads_index   // channel: [mandatory] [meta, bam, bai]
        disable_picard        // boolean: [optional]  [true]

    main:
        ch_versions = Channel.empty()
        ch_metrics  = Channel.empty()

        ch_meta_reads_index
        | multiMap { meta, reads, index ->
            meta_reads_index:   [meta, reads, index]
            meta_fasta:         [meta, WorkflowMain.getGenomeAttribute(meta.genome, 'fasta')]
            meta_fai:           [meta, WorkflowMain.getGenomeAttribute(meta.genome, 'fai')]
            meta_dict:          [meta, WorkflowMain.getGenomeAttribute(meta.genome, 'dict')]
        }
        | set {ch_to_qc}

        if (!disable_picard) {

            PICARD_COLLECTMULTIPLEMETRICS(
                ch_to_qc.meta_reads_index,
                ch_to_qc.meta_fasta,
                ch_to_qc.meta_fai
            )
            ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions.first())
            ch_metrics = ch_metrics.mix(PICARD_COLLECTMULTIPLEMETRICS.out.metrics)


            // ch_bam_bai_target_branched = ch_bam_bai_target.branch {
            //     hsmetrics  : it.size == 4 && it[3] != []
            //         return it
            //     wgsmetrics : true
            //         return [ it[0], it[1], it[2] ]
            // }

            // // WGS metrics
            // PICARD_COLLECTWGSMETRICS( ch_bam_bai_target_branched.wgsmetrics, ch_meta_fasta, ch_meta_fai, [] )
            // ch_versions = ch_versions.mix(PICARD_COLLECTWGSMETRICS.out.versions.first())
            // ch_metrics = ch_metrics.mix(PICARD_COLLECTWGSMETRICS.out.metrics)

            // // HS metrics
            // PICARD_COLLECTHSMETRICS( ch_bam_bai_target_branched.hsmetrics, ch_meta_fasta, ch_meta_fai, ch_fasta_dict)
            // ch_versions = ch_versions.mix(PICARD_COLLECTHSMETRICS.out.versions.first())
            // ch_metrics = ch_metrics.mix(PICARD_COLLECTHSMETRICS.out.metrics)

        }

        // SUBWORKFLOW: bam_stats_samtools
        // Run samtools QC modules
        // BAM_STATS_SAMTOOLS([meta, bam, bai])
        BAM_STATS_SAMTOOLS(
            ch_to_qc.meta_reads_index,
            ch_to_qc.meta_fasta
        )

        ch_metrics = ch_metrics.mix(
            BAM_STATS_SAMTOOLS.out.stats,
            BAM_STATS_SAMTOOLS.out.flagstat,
            BAM_STATS_SAMTOOLS.out.idxstats
        )
        ch_metrics.dump(tag: "BAM QC: metrics", pretty: true)
        ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)

    emit:
        metrics  = ch_metrics   // [[meta, metrics], [...], ...]
        versions = ch_versions  // [versions]
}
