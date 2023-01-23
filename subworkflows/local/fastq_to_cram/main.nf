#!/usr/bin/env nextflow

//
// Take fastq; align, postprocess and compress
//

// MODULES
include { ELPREP_SFM            } from "../../../modules/local/elprep/sfm/main.nf"
include { SAMTOOLS_INDEX        } from "../../../modules/nf-core/samtools/index/main.nf"
include { BIOBAMBAM_BAMSORMADUP } from "../../../modules/nf-core/biobambam/bamsormadup/main.nf"
// SUBWORKFLOWS
include { BAM_ARCHIVE       } from "../../local/bam_archive/main"
include { FASTQ_ALIGN_DNA   } from '../../nf-core/fastq_align_dna/main'
include { FASTA_INDEX_DNA   } from '../../nf-core/fasta_index_dna/main'



workflow FASTQ_TO_CRAM {
    take:
        ch_fastq            // channel: [mandatory] [meta, [fastq, ...]]
        ch_fasta_fai        // channel: [mandatory] [meta2, fasta, fai]
        aligner             // string:  [mandatory] aligner [bowtie2, bwamem, bwamem2, dragmap, snap]
        ch_aligner_index    // channel: [optional ] [meta2, aligner_index]
        postprocessor       // string:  [optional ] postprocessor [elprep, bamsormadup]

    main:

        ch_versions      = Channel.empty()
        ch_multiqc_files = Channel.empty()

        ch_fai        = ch_fasta_fai.map {meta, fasta, fai -> fai  }.collect()
        ch_fasta      = ch_fasta_fai.map {meta, fasta, fai -> fasta}.collect()

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // STEP: GENERATE ALIGNER INDEX
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        // Generate aligner index if not provided
        if ( ! ch_aligner_index ) {
            ch_aligner_index = FASTA_INDEX_DNA (
                ch_fasta_fai.map {meta, fasta, fai -> [meta, fasta]}, // channel: [meta, fasta]
                ch_fasta_fai.map {meta, fasta, fai -> [meta, []]},    // channel: [meta, altliftover] TODO: fix this
                aligner,                                              // string:  aligner
            )
            ch_versions = ch_versions.mix(FASTA_INDEX_DNA.out.versions)
        }
        ch_aligner_index.dump(tag: "FASTQ_TO_CRAM: aligner index",{FormattingService.prettyFormat(it)})

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // STEP: ALIGNMENT
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        ch_fastq.dump(tag: "FASTQ_TO_CRAM: reads to align",{FormattingService.prettyFormat(it)})

        // align fastq files per sample
        // ALIGNMENT([meta,fastq], index, sort)
        FASTQ_ALIGN_DNA(ch_fastq, ch_aligner_index, aligner, true)
        ch_versions = ch_versions.mix(FASTQ_ALIGN_DNA.out.versions)

        FASTQ_ALIGN_DNA.out.bam.dump(tag: "FASTQ_TO_CRAM: aligned bam", {FormattingService.prettyFormat(it)})

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // STEP: MARK DUPLICATES
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        // Gather bams per sample for merging and markdup
        ch_bam_per_sample = gather_split_files_per_sample(FASTQ_ALIGN_DNA.out.bam).dump(tag: "FASTQ_TO_CRAM: bam per sample",{FormattingService.prettyFormat(it)})

        ch_markdup_bam_bai = Channel.empty()
        switch (postprocessor) {
            case "elprep":
                // ELPREP_SFM([meta, bam]
                ELPREP_SFM(ch_bam_per_sample)

                // SAMTOOLS_INDEX([meta, bam])
                SAMTOOLS_INDEX(ELPREP_SFM.out.bam)

                ch_markdup_bam_bai = ELPREP_SFM.out.bam.join(SAMTOOLS_INDEX.out.bai)
                ch_multiqc_files = ch_multiqc_files.mix( ELPREP_SFM.out.metrics.map { meta, metrics -> return metrics} )
                ch_versions = ch_versions.mix(
                    ELPREP_SFM.out.versions,
                    SAMTOOLS_INDEX.out.versions
                )
                break
            case "bamsormadup":
                // BIOBAMBAM_BAMSORMADUP([meta, [bam, bam]], fasta)
                BIOBAMBAM_BAMSORMADUP(ch_bam_per_sample, ch_fasta)
                ch_markdup_bam_bai = BIOBAMBAM_BAMSORMADUP.out.bam.join(BIOBAMBAM_BAMSORMADUP.out.bam_index)
                ch_multiqc_files = ch_multiqc_files.mix( BIOBAMBAM_BAMSORMADUP.out.metrics.map { meta, metrics -> return metrics} )
                ch_versions = ch_versions.mix(BIOBAMBAM_BAMSORMADUP.out.versions)
                break
        }
        ch_markdup_bam_bai.dump(tag: "FASTQ_TO_CRAM: postprocessed bam", {FormattingService.prettyFormat(it)})

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // COMPRESSION AND CHECKSUM
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        BAM_ARCHIVE(
            ch_markdup_bam_bai,
            ch_fasta_fai,
        )
        ch_versions = ch_versions.mix(BAM_ARCHIVE.out.versions)

    emit:
        cram_crai       = BAM_ARCHIVE.out.cram_crai
        checksum        = BAM_ARCHIVE.out.checksum
        multiqc_files   = ch_multiqc_files
        versions        = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Gather split files per sample
def gather_split_files_per_sample(ch_files) {
    // Gather bam files per sample based on id
    ch_files.map {
        // set id to samplename, drop readgroup and count meta values
        meta, files ->
        gk = (meta.count ?: 1) * (meta.chunks ?: 1)
        return [
            groupKey(
                // replace id by samplename, drop readgroup meta, drop count and chunks
                meta - meta.subMap('id', 'readgroup', 'count', 'chunks') + [id: meta.samplename],
                gk
            ),
            files
        ]
    }
    .groupTuple(by:[0])
    .map { meta, files ->
        return [meta, files.flatten()]
    }
}
