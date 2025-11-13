#!/usr/bin/env nextflow

//
// FASTQ_ALIGN_RNA: Align RNAseq fastq files to a reference genome
//


include { STAR_ALIGN } from "../../../modules/nf-core/star/align/main.nf"

workflow FASTQ_ALIGN_RNA {
    take:
    ch_reads_aligner_index_gtf // channel: [mandatory] reads, aligner, index, gtf

    main:
    ch_bam = channel.empty()
    ch_reports = channel.empty()
    ch_versions = channel.empty()

    ch_reads_aligner_index_gtf
        .branch { meta, reads, aligner, index, gtf ->
            star: aligner == 'star'
            return [meta, reads, index, gtf]
            other: true
        }
        .set { ch_to_align }

    // Throw error for all samples with unsupported aligners
    ch_to_align.other.map { meta, _reads, aligner, _index, _fasta ->
        error("Unsupported aligner ${aligner} for sample ${meta.id}")
    }

    // Align fastq files to reference genome
    STAR_ALIGN(ch_to_align.star, "Illumina", "CMGG")
    // if aligner is STAR
    ch_bam = ch_bam.mix(STAR_ALIGN.out.bam)
    ch_reports = ch_reports.mix(
        STAR_ALIGN.out.log_final,
        STAR_ALIGN.out.log_progress,
        STAR_ALIGN.out.log_out,
    )
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())

    emit:
    bam      = ch_bam // channel: [ [meta], bam       ]
    reports  = ch_reports // channel: [ [meta], log       ]
    versions = ch_versions // channel: [ versions.yml      ]
}
