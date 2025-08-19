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
include { samplesheetToList } from 'plugin/nf-schema'


workflow CONSENSUS {
    take:
        ch_input_fastq                   // channel: [meta_with_readgroup, fastq] for SE/PE/duplex samples
        ch_genomes                       // map: reference genome files
        ch_umi_in_readname               // boolean
    main:

        ch_versions       = Channel.empty()

        ch_input_fastq
            .combine(ch_umi_in_readname)
            .map { meta, fq1, fq2, umi_flag ->
                tuple(meta, fq1, fq2, umi_flag)
            }
            .branch { _meta, _fq1, _fq2, umi_flag ->
                umi_in_readname: umi_flag == true
                umi_in_seq: umi_flag == false
            }
            .set  {ch_input_fastq_branch}


    // Readname branch

        // 1.1: FASTQ => uBAM
        ch_RN_uBAM = Channel.empty()
        RN_FQ  = ch_input_fastq_branch.umi_in_readname.map { meta, r1, r2, _f -> tuple(meta, [r1, r2]) }

        FASTQTOBAM_READNAME(RN_FQ)
        SAMTOOLS_INDEX(FASTQTOBAM_READNAME.out.bam)

        FASTQTOBAM_READNAME.out.bam
            .join(SAMTOOLS_INDEX.out.bai, by: 0)
            .map { meta, bam, bai -> tuple(meta, bam, bai) }
            .set { CH_FASTQTOBAM_WITH_BAI }

        FGBIO_COPYUMIFROMREADNAME(CH_FASTQTOBAM_WITH_BAI)
        ch_RN_uBAM = ch_RN_uBAM.mix(FGBIO_COPYUMIFROMREADNAME.out.bam)

        // 1.2: uBAM => Mapped BAM







/*
        // Seq branch => uBAM


        // 1.1: FASTQ => uBAM
        SEQ_FQ = ch_input_fastq_branch.umi_in_seq.map { meta, r1, r2, _f -> tuple(meta, [r1, r2]) }

        FASTQTOBAM_SEQ(SEQ_FQ)
*/

    emit:
        /*
        consensus_bam  = ch_consensus
        duplex_metrics = ch_duplex_metrics
        versions       = ch_versions
        */
        consensus_bam  = Channel.empty()
        versions       = Channel.empty()
}
