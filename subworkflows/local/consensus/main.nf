#!/usr/bin/env nextflow

include { FGBIO_COPYUMIFROMREADNAME               } from '../../../modules/nf-core/fgbio/copyumifromreadname/main'
include { FGBIO_CALLMOLECULARCONSENSUSREADS       } from '../../../modules/nf-core/fgbio/callmolecularconsensusreads/main'
include { FGBIO_FASTQTOBAM as FASTQTOBAM_READNAME } from '../../../modules/nf-core/fgbio/fastqtobam/main'
include { FGBIO_FASTQTOBAM as FASTQTOBAM_SEQ      } from '../../../modules/nf-core/fgbio/fastqtobam/main'
include { FGBIO_FILTERCONSENSUSREADS              } from '../../../modules/nf-core/fgbio/filterconsensusreads/main'
include { FGBIO_GROUPREADSBYUMI                   } from '../../../modules/nf-core/fgbio/groupreadsbyumi/main'
include { FGBIO_SORTBAM                           } from '../../../modules/nf-core/fgbio/sortbam/main'
include { FGBIO_ZIPPERBAMS                        } from '../../../modules/nf-core/fgbio/zipperbams/main'
include { SAMTOOLS_FASTQ                          } from '../../../modules/nf-core/samtools/fastq/main'
include { SAMTOOLS_INDEX                          } from '../../../modules/nf-core/samtools/index/main'
include { BWA_MEM                                 } from '../../../modules/nf-core/bwa/mem/main'


workflow CONSENSUS {
    take:
        ch_input_fastq                   // channel: [meta_with_readgroup, fastq] for SE/PE/duplex samples
        ch_genomes                       // map: reference genome files

    main:
        def ch_versions            = Channel.empty()
        def ch_fastq_readname      = Channel.empty()
        def ch_fastq_seq           = Channel.empty()
        def ch_ubam                = Channel.empty()
        def ch_ubam_with_bai       = Channel.empty()

    // 1.1: FASTQ => uBAM

        if (params.umi_in_readname) {

        // Case 1: UMI_in_readname


            ch_fastq_readname = ch_input_fastq.map { meta, r1, r2 -> tuple(meta, [r1, r2]) }

            FASTQTOBAM_READNAME(ch_fastq_readname)
            ch_versions = ch_versions.mix(FASTQTOBAM_READNAME.out.versions)

            SAMTOOLS_INDEX(FASTQTOBAM_READNAME.out.bam)
            ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

            FASTQTOBAM_READNAME.out.bam
                .join(SAMTOOLS_INDEX.out.bai, by: 0)
                .map { meta, bam, bai -> tuple(meta, bam, bai) }
                .set { ch_ubam_with_bai }

            FGBIO_COPYUMIFROMREADNAME(ch_ubam_with_bai)
            ch_versions = ch_versions.mix(FGBIO_COPYUMIFROMREADNAME.out.versions)

            ch_ubam = ch_ubam.mix(FGBIO_COPYUMIFROMREADNAME.out.bam)

        } else {

        // Case 2: UMI_in_sequence

            ch_fastq_seq = ch_input_fastq.map { meta, r1, r2 -> tuple(meta, [r1, r2]) }

            FASTQTOBAM_SEQ(ch_fastq_seq)
            ch_versions = ch_versions.mix(FASTQTOBAM_SEQ.out.versions)

            ch_ubam = ch_ubam.mix(FASTQTOBAM_SEQ.out.bam)
        }

    // 1.2: uBAM => Mapped BAM

        SAMTOOLS_FASTQ(ch_ubam)
        ch_versions = ch_versions.mix(SAMTOOLS_FASTQ.out.versions)

        def ch_bwa_index = Channel.value(file(params.genomes.GRCh38.bwamem))
        def ch_fasta     = Channel.value(file(params.genomes.GRCh38.fasta))

        BWA_MEM(
            SAMTOOLS_FASTQ.out.reads,
            ch_bwa_index,
            ch_fasta,
            false
        )
        ch_versions = ch_versions.mix(BWA_MEM.out.versions)

        FGBIO_ZIPPERBAMS(
            ch_ubam,
            BWA_MEM.out.bam,
            ch_fasta
        )
        ch_versions = ch_versions.mix(FGBIO_ZIPPERBAMS.out.versions)

    emit:
        /*
        consensus_bam  = ch_consensus
        duplex_metrics = ch_duplex_metrics
        versions       = ch_versions
        */
        consensus_bam  = FGBIO_ZIPPERBAMS.out.bam
        versions       = ch_versions
}
