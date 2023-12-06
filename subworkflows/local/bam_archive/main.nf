#!/usr/bin/env nextflow

include { SAMTOOLS_CONVERT } from "../../../modules/nf-core/samtools/convert/main"
include { MD5SUM           } from "../../../modules/nf-core/md5sum/main"

workflow BAM_ARCHIVE {
    take:
        ch_reads_index      // channel: [mandatory] [meta, bam/cram, bai,crai]
        ch_fasta_fai    // channel: [mandatory] [meta2, fasta, fai]

    main:
        ch_versions = Channel.empty()

        ch_meta_fai   = ch_fasta_fai.map {meta, fasta, fai -> [meta, fai]  }.collect()
        ch_meta_fasta = ch_fasta_fai.map {meta, fasta, fai -> [meta, fasta] }.collect()

        // separate bam and cram
        ch_reads_index = ch_reads_index.branch { meta, reads, index ->
            bam: reads.getExtension() == "bam"
                return [meta, reads, index]
            cram: reads.getExtension() == "cram"
                return [meta, reads, index]
        }

        // MODULE: samtools/convert
        // Compress bam to cram
        // SAMTOOLS CONVERT([meta, bam, bai], fasta, fai)
        SAMTOOLS_CONVERT(ch_reads_index.bam, ch_meta_fasta, ch_meta_fai)
        ch_versions = ch_versions.mix(SAMTOOLS_CONVERT.out.versions)


        ch_cram_crai = SAMTOOLS_CONVERT.out.alignment_index.mix(ch_reads_index.cram)

        // MODULE: MD5SUM
        // Generate md5sum for cram file
        // MD5SUM([meta, cram])
        MD5SUM(ch_cram_crai.map {meta, cram, crai -> [meta, cram]})
        ch_versions = ch_versions.mix(MD5SUM.out.versions)

    emit:
        cram_crai = ch_cram_crai            // [meta, cram, crai]
        checksum  = MD5SUM.out.checksum     // [meta, checksum]
        versions  = ch_versions             // [versions]
}
