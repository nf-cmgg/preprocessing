#!/usr/bin/env nextflow

//
// FASTQ_ALIGN_RNA: Align RNAseq fastq files to a reference genome
//


include { STAR_ALIGN                     } from "../../../modules/local/star/align/main"

workflow FASTQ_ALIGN_RNA {
    take:
        ch_reads_aligner_index_gtf      // channel: [mandatory] reads, aligner, index, gtf

    main:

        ch_bam_index    = Channel.empty()
        ch_bam          = Channel.empty()
        ch_reports      = Channel.empty()
        ch_versions     = Channel.empty()

        ch_reads_aligner_index_fasta.branch { meta, reads, aligner, index, fasta ->
            star : aligner == 'star'
                return [meta, reads, index, gtf]
            other   : true
        }
        .set{ch_to_align}

        // Throw error for all samples with unsupported aligners
        ch_to_align.other.map{ meta, reads, aligner, index, fasta ->
            error "Unsupported aligner ${aligner} for sample ${meta.id}"
        }

        // Align fastq files to reference genome
        STAR_ALIGN(ch_to_align.star) // if aligner is STAR
        ch_bam = ch_bam.mix(STAR_ALIGN.out.bam)
        ch_reports = ch_reports.mix(
            STAR_ALIGN.out.log_final,
            STAR_ALIGN.out.log_progress
            STAR_ALIGN.out.log_out
        )
        ch_versions = ch_versions.mix(STAR_ALIGN.out.versions)

    emit:
        bam         = ch_bam        // channel: [ [meta], bam       ]
        reports     = ch_reports    // channel: [ [meta], log       ]
        versions    = ch_versions   // channel: [ versions.yml      ]
}
