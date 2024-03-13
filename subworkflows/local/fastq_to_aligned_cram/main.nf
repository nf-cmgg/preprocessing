#!/usr/bin/env nextflow

//
// Take fastq; align, postprocess and compress
//

// MODULES
include { BIOBAMBAM_BAMSORMADUP } from "../../../modules/nf-core/biobambam/bamsormadup/main.nf"
include { SAMTOOLS_CONVERT      } from "../../../modules/nf-core/samtools/convert/main"
include { SAMTOOLS_SORMADUP     } from "../../../modules/nf-core/samtools/sormadup/main.nf"
include { SAMTOOLS_SORT         } from "../../../modules/nf-core/samtools/sort/main"

// SUBWORKFLOWS
include { FASTQ_ALIGN_DNA   } from '../../nf-core/fastq_align_dna/main'


workflow FASTQ_TO_CRAM {
    take:
        ch_meta_reads_alignerindex  // channel: [mandatory] [meta, [fastq, ...], aligner_index]
        aligner                     // string:  [mandatory] aligner [bowtie2, bwamem, bwamem2, dragmap, snap]
        markdup                     // string:  [optional ] markdup [bamsormadup, samtools, false]

    main:

        ch_versions      = Channel.empty()
        ch_multiqc_files = Channel.empty()

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // STEP: ALIGNMENT
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        ch_meta_reads_alignerindex.dump(tag: "FASTQ_TO_CRAM: reads to align",pretty: true)

        // align fastq files per sample
        // ALIGNMENT([meta,fastq], index, sort)
        FASTQ_ALIGN_DNA(
            ch_meta_reads_alignerindex,
            aligner,
            false
        )
        ch_versions = ch_versions.mix(FASTQ_ALIGN_DNA.out.versions)

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
        .groupTuple()
        .map {
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
        .groupTuple()
        .map { meta, files ->
            return [meta, files.flatten(), GenomeUtils.getGenomeAttribute(meta.genome, 'fasta')]
        }
        .set{ch_bam_fasta}
        ch_bam_fasta.dump(tag: "FASTQ_TO_CRAM: aligned bam per sample", pretty: true)

        ch_markdup_index = Channel.empty()

        switch (markdup) {
            case "bamsormadup":
                // BIOBAMBAM_BAMSORMADUP([meta, [bam, bam]], fasta)
                BIOBAMBAM_BAMSORMADUP(ch_bam_fasta)
                ch_markdup_index = ch_markdup_index.mix(BIOBAMBAM_BAMSORMADUP.out.bam.join(BIOBAMBAM_BAMSORMADUP.out.bam_index, failOnMismatch:true, failOnDuplicate:true))
                ch_multiqc_files = ch_multiqc_files.mix( BIOBAMBAM_BAMSORMADUP.out.metrics.map { meta, metrics -> return metrics} )
                ch_versions = ch_versions.mix(BIOBAMBAM_BAMSORMADUP.out.versions)
                break

            case "samtools":
                // SAMTOOLS_SORMADUP([meta, [bam, bam]], fasta)
                SAMTOOLS_SORMADUP(ch_bam_fasta)
                ch_markdup_index = ch_markdup_index.mix(SAMTOOLS_SORMADUP.out.cram.join(SAMTOOLS_SORMADUP.out.crai, failOnMismatch:true, failOnDuplicate:true))
                ch_multiqc_files = ch_multiqc_files.mix( SAMTOOLS_SORMADUP.out.metrics.map { meta, metrics -> return metrics} )
                ch_versions = ch_versions.mix(SAMTOOLS_SORMADUP.out.versions)
                break

            case ["false", false]:
                // Merge bam files and compress
                // SAMTOOLS_SORT([meta, [bam, bam], fasta])
                SAMTOOLS_SORT(ch_bam_fasta)
                ch_markdup_index = ch_markdup_index.mix(SAMTOOLS_SORT.out.cram.join(SAMTOOLS_SORT.out.crai, failOnMismatch:true, failOnDuplicate:true))
                ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)
                break
            default:
                error "markdup: ${markdup} not supported"

        }
        ch_markdup_index.dump(tag: "FASTQ_TO_CRAM: postprocessed bam", pretty: true)

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // COMPRESSION
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        ch_markdup_index
        .branch { meta, reads, index ->
            bam: reads.getExtension() == "bam"
                return [meta, reads, index]
            cram: reads.getExtension() == "cram"
                return [meta, reads, index]
        }
        .set {ch_markdup_index}

        ch_markdup_index.bam
        .map { meta, bam, bai ->
            bam_bai: [meta, bam, bai, GenomeUtils.getGenomeAttribute(meta.genome, 'fasta'), GenomeUtils.getGenomeAttribute(meta.genome, 'fai')]
        }
        .set {ch_bam_bai_fasta_fai}

        SAMTOOLS_CONVERT(ch_bam_bai_fasta_fai)

        ch_markdup_index.cram
        .mix(
            SAMTOOLS_CONVERT.out.cram.join(SAMTOOLS_CONVERT.out.crai, failOnMismatch:true, failOnDuplicate:true)
        )
        .set{ch_cram_crai}
        ch_cram_crai.dump(tag: "FASTQ_TO_CRAM: cram and crai", pretty: true)

    emit:
        cram_crai       = ch_cram_crai
        multiqc_files   = ch_multiqc_files
        versions        = ch_versions
}
