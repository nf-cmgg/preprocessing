#!/usr/bin/env nextflow

include { FGBIO_COPYUMIFROMREADNAME         } from '../../../modules/nf-core/fgbio/copyumifromreadname/main'
include { FGBIO_CALLMOLECULARCONSENSUSREADS } from '../../../modules/nf-core/fgbio/callmolecularconsensusreads/main'
include { FGBIO_FASTQTOBAM as FASTQTOBAM_READNAME           } from '../../../modules/nf-core/fgbio/fastqtobam/main'
include { FGBIO_FASTQTOBAM as FASTQTOBAM_SEQ           } from '../../../modules/nf-core/fgbio/fastqtobam/main'
include { FGBIO_FILTERCONSENSUSREADS        } from '../../../modules/nf-core/fgbio/filterconsensusreads/main'
include { FGBIO_GROUPREADSBYUMI             } from '../../../modules/nf-core/fgbio/groupreadsbyumi/main'
include { FGBIO_SORTBAM                     } from '../../../modules/nf-core/fgbio/sortbam/main'
include { FGBIO_ZIPPERBAMS                  } from '../../../modules/nf-core/fgbio/zipperbams/main'
include { SAMTOOLS_FASTQ                    } from '../../../modules/nf-core/samtools/fastq/main'
include { SAMTOOLS_INDEX                    } from '../../../modules/nf-core/samtools/index/main'
include { BWA_MEM                        } from '../../../modules/nf-core/bwa/mem/main'


workflow CONSENSUS {
    take:
        ch_input_fastq                   // channel: [meta_with_readgroup, fastq] for SE/PE/duplex samples
        ch_genomes                       // map: reference genome files

    main:
        def ch_versions       = Channel.empty()
        def readname_fastq     = Channel.empty()
        def seq_fastq = Channel.empty()
        def ch_fastq_to_uBAM   = Channel.empty()
        def ch_fastqtobam_with_bai = Channel.empty()

    // 1.1: FASTQ => uBAM

        if (params.umi_in_readname) {

        // Case 1: UMI_in_readname


            readname_fastq = ch_input_fastq.map { meta, r1, r2 -> tuple(meta, [r1, r2]) }

            FASTQTOBAM_READNAME(readname_fastq)
            ch_versions = ch_versions.mix(FASTQTOBAM_READNAME.out.versions)

            SAMTOOLS_INDEX(FASTQTOBAM_READNAME.out.bam)
            ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

            FASTQTOBAM_READNAME.out.bam
                .join(SAMTOOLS_INDEX.out.bai, by: 0)
                .map { meta, bam, bai -> tuple(meta, bam, bai) }
                .set { ch_fastqtobam_with_bai }

            FGBIO_COPYUMIFROMREADNAME(ch_fastqtobam_with_bai)
            ch_versions = ch_versions.mix(FGBIO_COPYUMIFROMREADNAME.out.versions)

            ch_fastq_to_uBAM = ch_fastq_to_uBAM.mix(FGBIO_COPYUMIFROMREADNAME.out.bam)

        } else {

        // Case 2: UMI_in_sequence

            seq_fastq = ch_input_fastq.map { meta, r1, r2 -> tuple(meta, [r1, r2]) }

            FASTQTOBAM_SEQ(seq_fastq)
            ch_versions = ch_versions.mix(FASTQTOBAM_SEQ.out.versions)

            ch_fastq_to_uBAM = ch_fastq_to_uBAM.mix(FASTQTOBAM_SEQ.out.bam)
        }

    // 1.2: uBAM => Mapped BAM
/*
    // Seq branch => uBAM


        // 1.1: FASTQ => uBAM
        SEQ_FQ = ch_input_fastq_branch.umi_in_seq.map { meta, r1, r2, _f -> tuple(meta, [r1, r2]) }

        FASTQTOBAM_SEQ(SEQ_FQ)

        // 1.2: uBAM -> Mapped BAM
*/

    emit:
        /*
        consensus_bam  = ch_consensus
        duplex_metrics = ch_duplex_metrics
        versions       = ch_versions
        */
        consensus_bam  = Channel.empty()
        versions       = ch_versions
}
