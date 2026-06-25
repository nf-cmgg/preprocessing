#!/usr/bin/env nextflow

//
// FASTQ_ALIGN_DNA: Align unaligned bam files to a reference genome
//

include { FGUMI_EXTRACT                   } from '../../../modules/nf-core/fgumi/extract/main.nf'
include { SNAPALIGNER_ALIGN as SNAP_ALIGN } from '../../../modules/nf-core/snapaligner/align/main'



workflow FASTQ_ALIGN_DNA {
    take:
    ch_reads_index // channel: [mandatory] reads, index

    main:

    ch_bam_index = channel.empty()
    ch_bam = channel.empty()
    ch_reports = channel.empty()


    ch_split_reads_index = ch_reads_index.multiMap { meta, reads, index ->
        reads: [meta, reads]
        index: [meta, index]
    }

    // Convert fastq to unaligned bam
    FGUMI_EXTRACT(
        ch_split_reads_index.reads.map { meta, reads ->
            def library = meta.library ?: 'default'
            return [meta, reads, library]
        }
    )

    FGUMI_EXTRACT.out.bam
        .join(ch_split_reads_index.index)
        .branch { meta, reads, index ->
            snap: meta.aligner == 'snap'
            return [meta, reads, index]
            other: true
        }
        .set { ch_to_align }

    // Throw error for all samples with unsupported aligners
    ch_to_align.other.map { meta, _reads, _index ->
        error("Unsupported aligner ${meta.aligner} for sample ${meta.id}")
    }

    // Align fastq files to reference genome and (optionally) sort
    // If aligner is snap
    SNAP_ALIGN(ch_to_align.snap)
    ch_bam = ch_bam.mix(SNAP_ALIGN.out.bam)
    ch_bam_index = ch_bam_index.mix(SNAP_ALIGN.out.bai)

    emit:
    bam       = ch_bam // channel: [ [meta], bam       ]
    bam_index = ch_bam_index // channel: [ [meta], csi/bai   ]
    reports   = ch_reports // channel: [ [meta], log       ]
}
