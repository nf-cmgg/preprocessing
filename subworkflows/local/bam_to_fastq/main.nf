#!/usr/bin/env nextflow

include { SAMTOOLS_FASTQ    } from '../../../modules/nf-core/modules/samtools/fastq/main'
include { SAMTOOLS_GETRG    } from '../../../modules/nf-core/modules/samtools/getrg/main'

workflow BAM_TO_FASTQ {
    take:
        ch_bam      // [meta, bam]
        ch_fasta    // fasta

    main:
        ch_versions = Channel.empty()
        ch_fastq  = Channel.empty()

        // Convert BAM/CRAM to FASTQ
        SAMTOOLS_FASTQ ( ch_bam )
        ch_fastq = SAMTOOLS_FASTQ.out.fastq
        ch_versions = ch_versions.mix(SAMTOOLS_FASTQ.out.versions)

        // Get RG info from BAM/CRAM
        SAMTOOLS_GETRG (ch_bam)
        ch_rg = SAMTOOLS_GETRG.out.readgroup.splitText().dump(tag: 'readgroups')
        ch_versions = ch_versions.mix(SAMTOOLS_GETRG.out.versions)

        // Add RG data to fastq meta

    emit:
        fastq    = ch_fastq    // [[meta, [fastq1, (fastq2)]]]
        versions = ch_versions // [versions]
}
