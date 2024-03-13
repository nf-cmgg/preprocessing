#!/usr/bin/env nextflow

//
// Take fastq; convert to ubam and compress

// MODULES
include { SAMTOOLS_CAT      } from '../../../modules/nf-core/samtools/cat/main'
include { SAMTOOLS_IMPORT   } from "../../../modules/nf-core/samtools/import/main"
include { MD5SUM            } from "../../../modules/nf-core/md5sum/main"

workflow FASTQ_TO_UCRAM {
    take:
        ch_fastq        // channel: [mandatory] [meta, [fastq, ...]]

    main:

        ch_versions = Channel.empty()

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // STEP: FASTQ TO BAM CONVERSION
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        ch_fastq
            .dump(tag: "FASTQ_TO_UCRAM: reads to convert",pretty: true)

        // SAMTOOLS_IMPORT([meta, fastq])
        SAMTOOLS_IMPORT(ch_fastq)
        ch_versions = ch_versions.mix(SAMTOOLS_IMPORT.first())

        SAMTOOLS_IMPORT.out.cram
        .map {
            // set id to samplename, drop readgroup and count meta values
            meta, files ->
            def gk = (meta.chunks as Integer ?: 1)
            return [
                groupKey(
                    // replace id by samplename, drop readgroup meta and chunks
                    meta - meta.subMap('id', 'readgroup', 'chunks') + [id: meta.samplename + ".unaligned"],
                    gk
                ),
                files
            ]
        }
        .groupTuple(by:[0])
        .dump(tag: "FASTQ_TO_UCRAM: unaligned cram per replicate",pretty: true)
        .map {
            meta, files ->
            def gk = (meta.count as Integer ?: 1)
            return [
                groupKey(
                    // drop count
                    meta - meta.subMap('count'),
                    gk
                ),
                files
            ]
        }
        .groupTuple(by:[0])
        .map { meta, files ->
            return [meta, files.flatten()]
        }
        .dump(tag: "FASTQ_TO_UCRAM: unaligned cram per sample",pretty: true)
        .set{ch_ubam_per_sample}

        // Merge bam files per sample
        SAMTOOLS_CAT(ch_ubam_per_sample)
        ch_versions = ch_versions.mix(SAMTOOLS_CAT.out.versions)

    emit:
        cram_crai = SAMTOOLS_CAT.out.cram    // [meta, cram]
        versions  = ch_versions              // versions

}
