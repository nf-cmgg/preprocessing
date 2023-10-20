#!/usr/bin/env nextflow

//
// Take fastq; align, postprocess and compress
//

// MODULES
include { SAMTOOLS_SORMADUP     } from "../../../modules/local/samtools/sormadup/main.nf"
include { BIOBAMBAM_BAMSORMADUP } from "../../../modules/nf-core/biobambam/bamsormadup/main.nf"
include { SAMTOOLS_SORTMERGE    } from "../../../modules/local/samtools/sortmerge/main"
// SUBWORKFLOWS
include { BAM_ARCHIVE       } from "../../local/bam_archive/main"
include { FASTQ_ALIGN_DNA   } from '../../nf-core/fastq_align_dna/main'


workflow FASTQ_TO_CRAM {
    take:
        ch_meta_reads   // channel: [mandatory] [meta, [fastq, ...]]
        aligner         // string:  [mandatory] aligner [bowtie2, bwamem, bwamem2, dragmap, snap]
        markdup         // string:  [optional ] markdup [bamsormadup, samtools, false]

    main:

        ch_versions      = Channel.empty()
        ch_multiqc_files = Channel.empty()

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // STEP: ALIGNMENT
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        ch_meta_reads.dump(tag: "FASTQ_TO_CRAM: reads to align",pretty: true)

        ch_meta_reads
        | multiMap { meta, reads ->
            reads: [meta, reads]
            index: [meta, WorkflowMain.getGenomeAttribute(meta.genome, aligner)]
        }
        | set {ch_to_align}

        // align fastq files per sample
        // ALIGNMENT([meta,fastq], index, sort)
        FASTQ_ALIGN_DNA(
            ch_to_align.reads,
            ch_to_align.index,
            aligner,
            false
        )
        ch_versions = ch_versions.mix(FASTQ_ALIGN_DNA.out.versions)

        FASTQ_ALIGN_DNA.out.bam.dump(tag: "FASTQ_TO_CRAM: aligned bam", pretty: true)

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // STEP: MARK DUPLICATES
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        FASTQ_ALIGN_DNA.out.bam.map {
            // set id to samplename, drop readgroup and count meta values
            meta, files ->
            def gk = (meta.chunks as Integer ?: 1)
            return [
                groupKey(
                    // replace id by samplename, drop readgroup meta and chunks
                    meta - meta.subMap('id', 'readgroup', 'chunks') + [id: meta.samplename],
                    gk
                ),
                files
            ]
        }
        | groupTuple(by:[0])
        | map {
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
        | groupTuple(by:[0])
        | map { meta, files ->
            return [meta, files.flatten()]
        }
        multiMap { meta, bam ->
            bam: [meta, bam]
            fasta: [meta, WorkflowMain.getGenomeAttribute(meta.genome, 'fasta')]
        }
        | set{ch_bam_per_sample}

        ch_markdup_bam_bai = Channel.empty()

        switch (markdup) {
            case "bamsormadup":
                // BIOBAMBAM_BAMSORMADUP([meta, [bam, bam]], fasta)
                BIOBAMBAM_BAMSORMADUP(ch_bam_per_sample.bam, ch_bam_per_sample.fasta)
                ch_markdup_bam_bai = BIOBAMBAM_BAMSORMADUP.out.bam.join(BIOBAMBAM_BAMSORMADUP.out.bam_index, failOnMismatch:true, failOnDuplicate:true)
                ch_multiqc_files = ch_multiqc_files.mix( BIOBAMBAM_BAMSORMADUP.out.metrics.map { meta, metrics -> return metrics} )
                ch_versions = ch_versions.mix(BIOBAMBAM_BAMSORMADUP.out.versions)
                break

            case "samtools":
                // SAMTOOLS_SORMADUP([meta, [bam, bam]], fasta)
                SAMTOOLS_SORMADUP(ch_bam_per_sample.bam, ch_bam_per_sample.fasta)
                ch_markdup_bam_bai = SAMTOOLS_SORMADUP.out.bam.join(SAMTOOLS_SORMADUP.out.bam_index, failOnMismatch:true, failOnDuplicate:true)
                ch_multiqc_files = ch_multiqc_files.mix( SAMTOOLS_SORMADUP.out.metrics.map { meta, metrics -> return metrics} )
                ch_versions = ch_versions.mix(SAMTOOLS_SORMADUP.out.versions)
                break

            case "false":
                // Merge bam files and compress
                // SAMTOOLS_MERGE([meta, [bam, bam]])
                SAMTOOLS_SORTMERGE(ch_bam_per_sample.bam)
                ch_markdup_bam_bai = SAMTOOLS_MERGE.out.bam.join(SAMTOOLS_SORTMERGE.out.bam_index, failOnMismatch:true, failOnDuplicate:true)
                ch_versions = ch_versions.mix(SAMTOOLS_SORTMERGE.out.versions)
                break

        }
        ch_markdup_bam_bai.dump(tag: "FASTQ_TO_CRAM: postprocessed bam", pretty: true)

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // COMPRESSION AND CHECKSUM
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        // BAM_ARCHIVE(
        //     ch_markdup_bam_bai,
        //     ch_fasta_fai,
        // )
        // ch_versions = ch_versions.mix(BAM_ARCHIVE.out.versions)

    emit:
        cram_crai       = BAM_ARCHIVE.out.cram_crai
        //checksum        = BAM_ARCHIVE.out.checksum
        multiqc_files   = ch_multiqc_files
        versions        = ch_versions
}
