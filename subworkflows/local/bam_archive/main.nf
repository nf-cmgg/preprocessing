#!/usr/bin/env nextflow

include { SAMTOOLS_CONVERT } from "../../../modules/nf-core/modules/samtools/convert/main"
include { MD5SUM           } from "../../../modules/nf-core/modules/md5sum/main"

workflow BAM_ARCHIVE {
    take:
        ch_bam_bai  // [meta, bam, bai]
        ch_fasta    // fasta
        ch_fai      // fasta.fai

    main:
        ch_versions = Channel.empty()

        // MODULE: samtools/convert
        // Compress bam to cram
        // SAMTOOLS CONVERT([meta, bam, bai], fasta, fai)
        SAMTOOLS_CONVERT(
            ch_bam_bai, ch_fasta, ch_fai
        )
        ch_versions = ch_versions.mix(SAMTOOLS_CONVERT.out.versions)

        // MODULE: MD5SUM
        // Generate md5sum for cram file
        // MD5SUM([meta, cram])
        MD5SUM(
            SAMTOOLS_CONVERT.out.alignment_index.map {
                meta, cram, crai -> return [meta, cram]
            }
        )
        ch_versions = ch_versions.mix(MD5SUM.out.versions)

    emit:
        cram_crai = SAMTOOLS_CONVERT.out.alignment_index  // [meta, cram, crai]
        checksums = MD5SUM.out.checksum                   // [meta, checksum]
        versions  = ch_versions                           // [versions]
}
