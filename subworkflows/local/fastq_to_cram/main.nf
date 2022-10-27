#!/usr/bin/env nextflow

//
// Take fastq; align, postprocess and compress
//

// MODULES
include { BIOBAMBAM_BAMSORMADUP } from "../../../modules/nf-core/biobambam/bamsormadup/main.nf"

// SUBWORKFLOWS
include { BAM_ARCHIVE       } from "../../local/bam_archive/main"
include { FASTQ_ALIGN_DNA   } from '../../nf-core/fastq_align_dna/main'
include { FASTA_INDEX_DNA   } from '../../nf-core/fasta_index_dna/main'



workflow FASTQ_TO_CRAM {
    take:
        ch_fastq_per_sample // channel: [mandatory] [meta, [fastq, ...]]
        ch_fasta_fai        // channel: [mandatory] [meta2, fasta, fai]
        aligner             // string:  [mandatory] aligner [bowtie2, bwamem, bwamem2, dragmap, snap]
        ch_aligner_index    // channel: [optional ] [meta2, aligner_index]

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
            ch_aligner_index = FASTA_INDEX_DNA ( ch_fasta_fai.map {meta, fasta, fai -> [meta, fasta]} )
            ch_versions = ch_versions.mix(FASTA_INDEX_DNA.out.versions)
        }
        ch_aligner_index.dump(tag: "FASTQ_TO_CRAM: aligner index",{FormattingService.prettyFormat(it)})

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // STEP: ALIGNMENT
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        ch_reads_to_map = Channel.empty()
        if (aligner == "snap" ) {
            ch_reads_to_map = ch_fastq_per_sample
        } else {
            ch_reads_to_map = ch_fastq_per_sample.map{meta, reads ->
                reads_files = meta.single_end ? reads : reads.sort().collate(2)
                return [meta, reads_files]
            }.transpose()
        }
        ch_reads_to_map.dump(tag: "FASTQ_TO_CRAM: reads to align",{FormattingService.prettyFormat(it)})

        // align fastq files per sample
        // ALIGNMENT([meta,fastq], index, sort)
        FASTQ_ALIGN_DNA(ch_reads_to_map, ch_aligner_index, aligner, true)
        ch_versions = ch_versions.mix(FASTQ_ALIGN_DNA.out.versions)

        FASTQ_ALIGN_DNA.out.bam.dump(tag: "FASTQ_TO_CRAM: aligned bam", {FormattingService.prettyFormat(it)})

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // STEP: MARK DUPLICATES
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        // Gather bams per sample for merging and markdup
        ch_markdup_bam_bai = Channel.empty()

        if (aligner == "snap") {
            ch_markdup_bam_bai = ch_markdup_bam_bai.mix(FASTQ_ALIGN_DNA.out.bam.join(FASTQ_ALIGN_DNA.out.bai))
        } else {
            ch_bam_per_sample = gather_split_files_per_sample(FASTQ_ALIGN_DNA.out.bam).dump(tag: "FASTQ_TO_CRAM: bam per sample",{FormattingService.prettyFormat(it)})

            // BIOBAMBAM_BAMSORMADUP([meta, [bam, bam]], fasta)
            BIOBAMBAM_BAMSORMADUP(ch_bam_per_sample, ch_fasta)
            ch_markdup_bam_bai = ch_markdup_bam_bai.mix(BIOBAMBAM_BAMSORMADUP.out.bam.join(BIOBAMBAM_BAMSORMADUP.out.bam_index))
            ch_multiqc_files = ch_multiqc_files.mix( BIOBAMBAM_BAMSORMADUP.out.metrics.map { meta, metrics -> return metrics} )
            ch_versions = ch_versions.mix(BIOBAMBAM_BAMSORMADUP.out.versions)
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
        // set id to filename without lane designation
        meta, files ->
        new_meta = meta.clone()
        new_meta.id = meta.samplename
        return [new_meta, files]
    }
    .groupTuple( by: [0])
    .map { meta, files ->
        return [meta, files.flatten()]
    }
}
