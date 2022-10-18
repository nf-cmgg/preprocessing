#!/usr/bin/env nextflow

include { SAMTOOLS_FASTQ    } from '../../../modules/nf-core/samtools/fastq/main'
include { SAMTOOLS_GETRG    } from '../../../modules/nf-core/samtools/getrg/main'

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

        // Get RG info from BAM/CRAM and parse into maps
        SAMTOOLS_GETRG (ch_bam)
        ch_rg = SAMTOOLS_GETRG.out.readgroup
            .splitText()
            .map{ meta, rg ->
                meta.readgroup = rg.stripTrailing().split('\t').tail().collectEntries{ [it.split(':')[0], it.split(':')[1]] }
                meta.readgroup['SM'] = meta.samplename
                return [meta]
            }
        ch_versions = ch_versions.mix(SAMTOOLS_GETRG.out.versions)

        // Add RG data to fastq meta
        // thanks to @Midnighter for the utility function
        ch_fastq_with_rg = CustomChannelOperators.joinOnKeys(ch_rg, ch_fastq, "samplename").dump(tag: 'fastq with RG', {FormattingService.prettyFormat(it)})
    emit:
        fastq    = ch_fastq_with_rg    // [[meta, [fastq1, (fastq2)]]]
        versions = ch_versions // [versions]
}
