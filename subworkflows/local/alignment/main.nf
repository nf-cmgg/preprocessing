#!/usr/bin/env nextflow

//
// ALIGNMENT
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { BIOBAMBAM_BAMSORMADUP             } from "../../../modules/nf-core/modules/biobambam/bamsormadup/main"
include { BOWTIE2_ALIGN                     } from "../../../modules/nf-core/modules/bowtie2/align/main"
include { BWA_MEM as BWAMEM1_MEM            } from '../../../modules/nf-core/modules/bwa/mem/main'
include { BWAMEM2_MEM                       } from '../../../modules/nf-core/modules/bwamem2/mem/main'
include { DRAGMAP_ALIGN                     } from '../../../modules/nf-core/modules/dragmap/align/main'
include { SNAPALIGNER_ALIGN as SNAP_ALIGN   } from '../../../modules/nf-core/modules/snapaligner/align/main'


workflow ALIGNMENT {
    take:
        ch_reads            // channel: [mandatory] meta, reads
        ch_aligner_index    // channel: [mandatory] aligner index
        sort                // boolean: [mandatory] true -> sort, false -> don't sort

    main:

        ch_versions = Channel.empty()

        // Align fastq files to reference genome and (optionally) sort
        BOWTIE2(ch_reads, ch_aligner_index, false, sort)// if aligner is bowtie2
        BWAMEM1_MEM(ch_reads,   ch_map_index, sort)     // If aligner is bwa-mem
        BWAMEM2_MEM(ch_reads,   ch_map_index, sort)     // If aligner is bwa-mem2
        DRAGMAP_ALIGN(ch_reads, ch_map_index, sort)     // If aligner is dragmap
        SNAP_ALIGN(ch_reads, ch_map_index)              // If aligner is snap

        ch_versions = ch_versions.mix(
            BOWTIE2.out.versions,
            BWAMEM1_MEM.versions,
            BWAMEM2_MEM.versions,
            DRAGMAP_ALIGN.versions,
            SNAP_ALIGN.versions
        )

        // Get the bam files from the aligner
        // Only one aligner is run
        ch_bam = Channel.empty()
        ch_bam = ch_bam.mix(
            BOWTIE2.out.bam,
            BWAMEM1_MEM.out.bam,
            BWAMEM2_MEM.out.bam,
            DRAGMAP_ALIGN.out.bam,
            SNAP_ALIGN.out.bam,
        ).groupTuple()

        // MODULE: bamsormadup
        // Merge, sort and mark duplicates
        // BIOBAMBAM_BAMSORMADUP([meta, bam1, bam2, ...], fasta)
        BIOBAMBAM_BAMSORMADUP(ch_bam, [])
        ch_markdup_bam_bai = BIOBAMBAM_BAMSORMADUP.out.bam.join(BIOBAMBAM_BAMSORMADUP.out.bam_index)
        ch_markdup_metrics = BIOBAMBAM_BAMSORMADUP.out.metrics
        ch_versions = ch_versions.mix(BIOBAMBAM_BAMSORMADUP.out.versions)

    emit:
        bam_bai         = ch_markdup_bam_bai // channel: [ [meta], bam, bai ]
        markdup_metrics = ch_markdup_metrics // channel: [ [meta], metrics ]
        versions        = ch_versions        // channel: [ versions.yml ]
}

