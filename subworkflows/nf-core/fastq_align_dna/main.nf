#!/usr/bin/env nextflow

//
// FASTQ_ALIGN_DNA: Align fastq files to a reference genome
//


include { BOWTIE2_ALIGN                     } from "../../../modules/nf-core/bowtie2/align/main"
include { BWA_MEM as BWAMEM1_MEM            } from '../../../modules/nf-core/bwa/mem/main'
include { BWAMEM2_MEM as BWAMEM2_MEM        } from '../../../modules/nf-core/bwamem2/mem/main'
include { DRAGMAP_ALIGN                     } from "../../../modules/nf-core/dragmap/align/main"
include { SNAPALIGNER_ALIGN as SNAP_ALIGN   } from '../../../modules/nf-core/snapaligner/align/main'



workflow FASTQ_ALIGN_DNA {
    take:
        ch_reads_aligner_index_fasta    // channel: [mandatory] reads, aligner, index, fasta
        sort                            // boolean: [mandatory] true -> sort, false -> don't sort

    main:

        ch_bam_index    = Channel.empty()
        ch_bam          = Channel.empty()
        ch_reports      = Channel.empty()
        ch_versions     = Channel.empty()

        ch_reads_aligner_index_fasta.branch { meta, reads, aligner, index, fasta ->
            bowtie2 : aligner == 'bowtie2'
                return [meta, reads, index, fasta]
            bwamem  : aligner == 'bwamem'
                return [meta, reads, index, fasta]
            bwamem2 : aligner == 'bwamem2'
                return [meta, reads, index, fasta]
            dragmap : aligner == 'dragmap'
                return [meta, reads, index, fasta]
            snap    : aligner == 'snap'
                return [meta, reads, index]
            other   : true
        }
        .set{ch_to_align}

        // Throw error for all samples with unsupported aligners
        ch_to_align.other.map{ meta, reads, aligner, index, fasta ->
            error "Unsupported aligner ${aligner} for sample ${meta.id}"
        }

        // Align fastq files to reference genome and (optionally) sort
        BOWTIE2_ALIGN(ch_to_align.bowtie2, false, sort) // if aligner is bowtie2
        ch_bam = ch_bam.mix(BOWTIE2_ALIGN.out.bam)
        ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions)

        BWAMEM1_MEM  (ch_to_align.bwamem, sort)         // If aligner is bwa-mem
        ch_bam = ch_bam.mix(BWAMEM1_MEM.out.bam)
        ch_bam_index = ch_bam_index.mix(BWAMEM1_MEM.out.csi)
        ch_versions = ch_versions.mix(BWAMEM1_MEM.out.versions)

        BWAMEM2_MEM  (ch_to_align.bwamem2, sort)        // If aligner is bwa-mem2
        ch_bam = ch_bam.mix(BWAMEM2_MEM.out.bam)
        ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions)

        DRAGMAP_ALIGN(ch_to_align.dragmap, sort)        // If aligner is dragmap
        ch_bam = ch_bam.mix(DRAGMAP_ALIGN.out.bam)
        ch_reports = ch_reports.mix(DRAGMAP_ALIGN.out.log)
        ch_versions = ch_versions.mix(DRAGMAP_ALIGN.out.versions)

        SNAP_ALIGN(ch_to_align.snap)                    // If aligner is snap
        ch_bam = ch_bam.mix(SNAP_ALIGN.out.bam)
        ch_bam_index.mix(SNAP_ALIGN.out.bai)
        ch_versions = ch_versions.mix(SNAP_ALIGN.out.versions)

    emit:
        bam         = ch_bam        // channel: [ [meta], bam       ]
        bam_index   = ch_bam_index  // channel: [ [meta], csi/bai   ]
        reports     = ch_reports    // channel: [ [meta], log       ]
        versions    = ch_versions   // channel: [ versions.yml      ]
}
