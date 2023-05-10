#!/usr/bin/env nextflow

//
// Take fastq; convert to ubam and compress

// MODULES
include { SAMTOOLS_CAT      } from '../../../modules/nf-core/samtools/cat/main'
include { FGBIO_FASTQTOBAM  } from "../../../modules/nf-core/fgbio/fastqtobam/main"

// SUBWORKFLOWS
include { BAM_ARCHIVE       } from "../../local/bam_archive/main"

workflow FASTQ_TO_UCRAM {
    take:
        ch_fastq        // channel: [mandatory] [meta, [fastq, ...]]
        ch_fasta_fai    // channel: [mandatory] [meta2, fasta, fai]

    main:

        ch_versions = Channel.empty()

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // STEP: FASTQ TO BAM CONVERSION
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        ch_fastq
            .dump(tag: "FASTQ_TO_UCRAM: reads to convert",{FormattingService.prettyFormat(it)})

        // FGBIO_FASTQTOBAM([meta, fastq])
        FGBIO_FASTQTOBAM(ch_fastq)
        ch_versions = ch_versions.mix(FGBIO_FASTQTOBAM.out.versions)

        FGBIO_FASTQTOBAM.out.bam.map {
            // set id to samplename, drop readgroup and count meta values
            meta, files ->
            def gk = (meta.chunks ?: 1)
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
        .dump(tag: "FASTQ_TO_UCRAM: unaligned bam per replicate",{FormattingService.prettyFormat(it)})
        .map {
            meta, files ->
            def gk = (meta.count ?: 1)
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
        .dump(tag: "FASTQ_TO_UCRAM: unaligned bam per sample",{FormattingService.prettyFormat(it)})
        .set{ch_ubam_per_sample}

        // Merge bam files per sample
        SAMTOOLS_CAT(ch_ubam_per_sample)

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // COMPRESSION AND CHECKSUM
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        BAM_ARCHIVE(SAMTOOLS_CAT.out.bam.map { meta,bam -> [meta, bam, []]}, ch_fasta_fai)
        ch_versions = ch_versions.mix(BAM_ARCHIVE.out.versions)

    emit:
        checksum  = BAM_ARCHIVE.out.checksum    // [meta, checksum]
        cram_crai = BAM_ARCHIVE.out.cram_crai   // [meta, cram, crai]
        versions  = ch_versions                 // versions

}
