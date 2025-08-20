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
include { BWAMEM2_MEM                             } from '../../../modules/nf-core/bwamem2/mem/main'
include { BOWTIE2_ALIGN                           } from '../../../modules/nf-core/bowtie2/align/main'


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

        def ch_reads_aligner_index_fasta = SAMTOOLS_FASTQ.out.reads.map { meta, reads ->
            def gd    = (meta.genome_data instanceof Map) ? meta.genome_data : [:]
            def alg   = (meta.aligner ?: 'bwamem')
            def fasta = file(gd.fasta, checkIfExists: true)
            def index = file(gd[alg],  checkIfExists: true)
            tuple(meta, reads, alg, index, fasta)
        }

        ch_reads_aligner_index_fasta.branch { meta, reads, alg, index, fasta ->
            bwamem  : alg == 'bwamem'  ; return [meta, reads, index, fasta]
            bwamem2 : alg == 'bwamem2' ; return [meta, reads, index, fasta]
            bowtie2 : alg == 'bowtie2' ; return [meta, reads, index, fasta]
            other   : true
        }.set { ch_to_map }

        BWA_MEM(ch_to_map.bwamem, false)
        BWAMEM2_MEM(ch_to_map.bwamem2, false)
        BOWTIE2_ALIGN(ch_to_map.bowtie2, false, false)

        ch_versions = ch_versions.mix(BWA_MEM.out.versions)
        ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions)
        ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions)

        def ch_mapped_bam = Channel.empty()
        ch_mapped_bam = ch_mapped_bam.mix(BWA_MEM.out.bam)
        ch_mapped_bam = ch_mapped_bam.mix(BWAMEM2_MEM.out.bam)
        ch_mapped_bam = ch_mapped_bam.mix(BOWTIE2_ALIGN.out.bam)

        def ch_fasta_by_meta = ch_reads_aligner_index_fasta.map { meta, _r, _a, _i, fasta -> tuple(meta, fasta) }

        FGBIO_ZIPPERBAMS(
            ch_ubam,
            ch_mapped_bam,
            ch_fasta_by_meta
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
