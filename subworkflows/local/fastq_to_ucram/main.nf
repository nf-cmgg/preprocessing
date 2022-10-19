#!/usr/bin/env nextflow

//
// Take fastq; convert to ubam and compress

// MODULES
include { CAT_FASTQ            } from "../../../modules/nf-core/cat/fastq/main"
include { BIOBAMBAM_BAMSORMADUP       } from "../../../modules/nf-core/biobambam/bamsormadup/main"
include { FGBIO_FASTQTOBAM            } from "../../../modules/nf-core/fgbio/fastqtobam/main"

// SUBWORKFLOWS
include { BAM_ARCHIVE       } from "../../local/bam_archive/main"

workflow FASTQ_TO_UCRAM {
    take:
        ch_fastq_per_sample // channel: [mandatory] [meta, [fastq, ...]]
        ch_fasta_fai        // channel: [mandatory] [meta2, fasta, fai]

    main:

        ch_versions = Channel.empty()

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // STEP: FASTQ TO BAM CONVERSION
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        // Convert non-standard fastq data (e.g. non-human, non-DNA, ...) to BAM
        // CAT_FASTQ([meta, fastq])
        // Merge split fastqs to simplify converting to uBAM
        CAT_FASTQ(ch_fastq_per_sample)
        ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)

        // FGBIO_FASTQTOBAM([meta, fastq])
        FGBIO_FASTQTOBAM(CAT_FASTQ.out.reads)
        ch_versions = ch_versions.mix(FGBIO_FASTQTOBAM.out.versions)

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // COMPRESSION AND CHECKSUM
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        BAM_ARCHIVE(FGBIO_FASTQTOBAM.out.bam, ch_fasta_fai)
        ch_versions = ch_versions.mix(BAM_ARCHIVE.out.versions)

    emit:
        checksum  = BAM_ARCHIVE.out.checksum    // [meta, checksum]
        cram_crai = BAM_ARCHIVE.out.cram_crai   // [meta, cram, crai]
        versions  = ch_versions                 // versions

}
