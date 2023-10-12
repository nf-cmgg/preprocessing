#!/usr/bin/env nextflow

include { SAMTOOLS_COLLATEFASTQ } from '../../../modules/nf-core/samtools/collatefastq/main'
include { SAMTOOLS_GETRG        } from '../../../modules/nf-core/samtools/getrg/main'

workflow BAM_TO_FASTQ {
    take:
        ch_bam       // channel: [mandatory] [meta, bam]
        ch_fasta_fai // channel: [mandatory] [meta2, fasta, fai]

    main:
        ch_versions = Channel.empty()
        ch_fastq  = Channel.empty()

        ch_fai        = ch_fasta_fai.map {meta, fasta, fai -> fai          }.collect()
        ch_fasta      = ch_fasta_fai.map {meta, fasta, fai -> fasta        }.collect()
        ch_meta_fasta = ch_fasta_fai.map {meta, fasta, fai -> [meta, fasta]}.collect()

        // Convert BAM/CRAM to FASTQ
        SAMTOOLS_COLLATEFASTQ ( ch_bam, ch_meta_fasta, false )

        // Only keep the R1 and R2 FASTQ files
        ch_fastq = SAMTOOLS_COLLATEFASTQ.out.fastq
        ch_versions = ch_versions.mix(SAMTOOLS_COLLATEFASTQ.out.versions)

        // Get RG info from BAM/CRAM and parse into maps
        SAMTOOLS_GETRG (ch_bam)
        ch_rg = SAMTOOLS_GETRG.out.readgroup
            .splitText()
            .map{ meta, rg ->
                meta.readgroup = rg.stripTrailing().split('\t').tail().collectEntries{ [it.split(':')[0], it.split(':')[1]] }
                meta.readgroup['SM'] = meta.samplename ?: meta.id
                return [meta]
            }
        ch_versions = ch_versions.mix(SAMTOOLS_GETRG.out.versions)

        // Add RG data to fastq meta
        // thanks to @Midnighter for the utility function
        ch_fastq_with_rg = CustomChannelOperators.joinOnKeys(ch_rg, ch_fastq, "samplename").dump(tag: 'fastq with RG', pretty: true)
    emit:
        fastq    = ch_fastq_with_rg // [meta, [fastq1, (fastq2)]]
        versions = ch_versions      // [versions]
}
