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
        ch_reads_aligner_index_fasta    // channel: [mandatory] reads, aligner, aligner_index, fasta file
        sort                            // boolean: [mandatory] true -> sort, false -> don't sort

    main:

        ch_bam_index    = Channel.empty()
        ch_bam          = Channel.empty()
        ch_reports      = Channel.empty()
        ch_versions     = Channel.empty()

        ch_reads_aligner_index_fasta.branch {
            bowtie2 : it.filter{ it.meta.aligner == 'bowtie2' }
            bwamem  : it.filter{ it.meta.aligner == 'bwamem' }
            bwamem2 : it.filter{ it.meta.aligner == 'bwamem2' }
            dragmap : it.filter{ it.meta.aligner == 'dragmap' }
            snap    : it.filter{ it.meta.aligner == 'snap' }
            other      : true
        }
        .set{ch_to_align}

        // Align fastq files to reference genome and (optionally) sort
        BOWTIE2_ALIGN(ch_to_align.bowtie2, false, sort) // if aligner is bowtie2
        ch_bam = ch_bam.mix(BOWTIE2_ALIGN.out.bam)
        ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions)

        BWAMEM1_MEM  (ch_to_align.bwamem, sort)        // If aligner is bwa-mem
        ch_bam = ch_bam.mix(BWAMEM1_MEM.out.bam)
        ch_bam_index = ch_bam_index.mix(BWAMEM1_MEM.out.csi)
        ch_versions = ch_versions.mix(BWAMEM1_MEM.out.versions)

        BWAMEM2_MEM  (ch_to_align.bwamem2, sort)       // If aligner is bwa-mem2
        ch_bam = ch_bam.mix(BWAMEM2_MEM.out.bam)
        ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions)

        DRAGMAP_ALIGN(ch_to_align.dragmap, sort)       // If aligner is dragmap
        ch_bam = ch_bam.mix(DRAGMAP_ALIGN.out.bam)
        ch_reports = ch_reports.mix(DRAGMAP_ALIGN.out.log)
        ch_versions = ch_versions.mix(DRAGMAP_ALIGN.out.versions)

        SNAP_ALIGN(ch_to_align.snap.map{ meta, reads, index, fasta -> return [meta, reads, index]}) // If aligner is snap
        ch_bam = ch_bam.mix(SNAP_ALIGN.out.bam)
        ch_bam_index.mix(SNAP_ALIGN.out.bai)
        ch_versions = ch_versions.mix(SNAP_ALIGN.out.versions)

    emit:
        bam         = ch_bam        // channel: [ [meta], bam       ]
        bam_index   = ch_bam_index  // channel: [ [meta], csi/bai   ]
        reports     = ch_reports    // channel: [ [meta], log       ]
        versions    = ch_versions   // channel: [ versions.yml      ]
}
