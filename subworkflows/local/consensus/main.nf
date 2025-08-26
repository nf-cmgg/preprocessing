#!/usr/bin/env nextflow

include { FGBIO_COPYUMIFROMREADNAME                  } from '../../../modules/nf-core/fgbio/copyumifromreadname/main'
include { FGBIO_CALLMOLECULARCONSENSUSREADS          } from '../../../modules/nf-core/fgbio/callmolecularconsensusreads/main'
include { FGBIO_FASTQTOBAM as FASTQTOBAM_READNAME    } from '../../../modules/nf-core/fgbio/fastqtobam/main'
include { FGBIO_FASTQTOBAM as FASTQTOBAM_SEQ         } from '../../../modules/nf-core/fgbio/fastqtobam/main'
include { FGBIO_FILTERCONSENSUSREADS                 } from '../../../modules/nf-core/fgbio/filterconsensusreads/main'
include { FGBIO_GROUPREADSBYUMI                      } from '../../../modules/nf-core/fgbio/groupreadsbyumi/main'
include { FGBIO_SORTBAM                              } from '../../../modules/nf-core/fgbio/sortbam/main'
include { FGBIO_ZIPPERBAMS as FGBIO_ZIPPERBAMS_RAW   } from '../../../modules/nf-core/fgbio/zipperbams/main'
include { FGBIO_ZIPPERBAMS as FGBIO_ZIPPERBAMS_CONS  } from '../../../modules/nf-core/fgbio/zipperbams/main'
include { SAMTOOLS_FASTQ as SAMTOOLS_FASTQ_RAW       } from '../../../modules/nf-core/samtools/fastq/main'
include { SAMTOOLS_FASTQ as SAMTOOLS_FASTQ_CONS      } from '../../../modules/nf-core/samtools/fastq/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_UMI       } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_CONSENSUS } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_CONVERT as SAMTOOLS_CONVERT       } from '../../../modules/nf-core/samtools/convert/main'
include { FASTQ_ALIGN_DNA as FASTQ_ALIGN_DNA_RAW     } from '../../../subworkflows/nf-core/fastq_align_dna/main'
include { FASTQ_ALIGN_DNA as FASTQ_ALIGN_DNA_CONS    } from '../../../subworkflows/nf-core/fastq_align_dna/main'

workflow CONSENSUS {
    take:
        ch_umi_fastq                   // channel: [meta, fastq1, fastq2]

    main:
        def ch_versions            = Channel.empty()
        def ch_ubam                = Channel.empty()

    // 1.1: FASTQ => uBAM

        def ch_fastq = ch_umi_fastq
            .map { meta, r1, r2 -> tuple(meta, [r1, r2]) }
            .branch { meta, _fastq ->
                readname: meta['umi_type'] == 'readname'
                seq:      meta['umi_type'] == 'seq'
            }

        // Case 1: UMI_in_readname
        if (ch_fastq.readname) {
            FASTQTOBAM_READNAME(ch_fastq.readname)
            ch_versions = ch_versions.mix(FASTQTOBAM_READNAME.out.versions)

            SAMTOOLS_INDEX_UMI(FASTQTOBAM_READNAME.out.bam)
            ch_versions = ch_versions.mix(SAMTOOLS_INDEX_UMI.out.versions)

            FASTQTOBAM_READNAME.out.bam
                .join(SAMTOOLS_INDEX_UMI.out.bai, by: 0)
                .map { meta, bam, bai -> tuple(meta, bam, bai) }
                .set { ch_ubam_with_bai }

            FGBIO_COPYUMIFROMREADNAME(ch_ubam_with_bai)
            ch_versions = ch_versions.mix(FGBIO_COPYUMIFROMREADNAME.out.versions)

            ch_ubam = ch_ubam.mix(FGBIO_COPYUMIFROMREADNAME.out.bam)
        }

        // Case 2: UMI_in_sequence

        if (ch_fastq.seq) {
            FASTQTOBAM_SEQ(ch_fastq.seq)
            ch_versions = ch_versions.mix(FASTQTOBAM_SEQ.out.versions)

            ch_ubam = ch_ubam.mix(FASTQTOBAM_SEQ.out.bam)
        }

    // 1.2: uBAM => Mapped BAM

        SAMTOOLS_FASTQ_RAW(ch_ubam, true)

        ch_versions = ch_versions.mix(SAMTOOLS_FASTQ_RAW.out.versions)

        def ch_reads_aligner_index_fasta = SAMTOOLS_FASTQ_RAW.out.interleaved.map { meta, reads ->
            def gd    = (meta.genome_data instanceof Map) ? meta.genome_data : [:]
            def alg   = (meta.aligner ?: 'bwamem')
            def fasta = file(gd.fasta, checkIfExists: true)
            def index = file(gd[alg],  checkIfExists: true)
            tuple(meta, reads, alg, index, fasta)
        }

        FASTQ_ALIGN_DNA_RAW(ch_reads_aligner_index_fasta, false)
        ch_versions = ch_versions.mix(FASTQ_ALIGN_DNA_RAW.out.versions)

        def ch_mapped_bam = FASTQ_ALIGN_DNA_RAW.out.bam

        def ch_fasta_by_meta = ch_reads_aligner_index_fasta.map { meta, _r, _a, _i, fasta -> tuple(meta, fasta) }

        def ch_dict_by_meta = ch_reads_aligner_index_fasta.map { meta, _r, _a, _i, _fasta ->
            def dict = file(meta.genome_data.dict, checkIfExists: true)
            tuple(meta, dict)
        }
        ch_ubam
            .join(ch_mapped_bam, by:0)
            .join(ch_fasta_by_meta, by:0)
            .join(ch_dict_by_meta, by:0)
            .map { meta, ubam, mapped_bam, fasta, dict -> tuple(meta, ubam, mapped_bam, fasta, dict) }
            .set { ch_zipperbam }

        FGBIO_ZIPPERBAMS_RAW(ch_zipperbam)

        ch_versions = ch_versions.mix(FGBIO_ZIPPERBAMS_RAW.out.versions)

    // 1.3: Mapped BAM => Grouped BAM


        def ch_strategy = Channel.value(params.umi_group_strategy)

        FGBIO_GROUPREADSBYUMI(
            FGBIO_ZIPPERBAMS_RAW.out.bam,
            ch_strategy
        )

        ch_versions     = ch_versions.mix(FGBIO_GROUPREADSBYUMI.out.versions)
        def ch_grouped_bam = FGBIO_GROUPREADSBYUMI.out.bam

    // 2(b).1: GroupedBam -> Filtered Consensus uBam
        def call_min_reads = Channel.value(params.callmolecularconsensusreads_min_reads)
        def call_min_baseq = Channel.value(params.callmolecularconsensusreads_min_baseq)

        FGBIO_CALLMOLECULARCONSENSUSREADS(
            ch_grouped_bam,
            call_min_reads,
            call_min_baseq
        )

        ch_versions = ch_versions.mix(FGBIO_CALLMOLECULARCONSENSUSREADS.out.versions)

        def ch_input_filterconsensusreads = FGBIO_CALLMOLECULARCONSENSUSREADS.out.bam.map {meta, bam ->
            def fasta = file(meta.genome_data.fasta, checkIfExists: true)
            tuple(meta, bam, fasta)
        }

        def filter_min_reads = Channel.value(params.filterconsensusreads_min_reads)
        def filter_min_baseq = Channel.value(params.filterconsensusreads_min_baseq)
        def filter_min_base_error_rate = Channel.value(params.filterconsensusreads_min_base_error_rate)

        FGBIO_FILTERCONSENSUSREADS(
            ch_input_filterconsensusreads,
            filter_min_reads,
            filter_min_baseq,
            filter_min_base_error_rate
        )

        ch_versions = ch_versions.mix(FGBIO_FILTERCONSENSUSREADS.out.versions)

        ch_filtered_uBam = FGBIO_FILTERCONSENSUSREADS.out.bam

    // 2(b).2: Consensus Filtered uBam -> Consensus Mapped & Filtered BAM

        SAMTOOLS_FASTQ_CONS(ch_filtered_uBam, true)
        ch_versions = ch_versions.mix(SAMTOOLS_FASTQ_CONS.out.versions)

        def ch_cons_reads_aligner_index_fasta = SAMTOOLS_FASTQ_CONS.out.interleaved.map { meta, reads ->
            def gd    = (meta.genome_data instanceof Map) ? meta.genome_data : [:]
            def alg   = (meta.aligner ?: 'bwamem')
            def fasta = file(gd.fasta, checkIfExists: true)
            def index = file(gd[alg],  checkIfExists: true)
            tuple(meta, reads, alg, index, fasta)
        }

        FASTQ_ALIGN_DNA_CONS(ch_cons_reads_aligner_index_fasta, false)
        ch_versions = ch_versions.mix(FASTQ_ALIGN_DNA_CONS.out.versions)

        def ch_cons_mapped_bam = FASTQ_ALIGN_DNA_CONS.out.bam
        def ch_cons_fasta_by_meta = ch_cons_reads_aligner_index_fasta.map { meta, _r, _a, _i, fasta -> tuple(meta, fasta) }
        def ch_cons_dict_by_meta  = ch_cons_reads_aligner_index_fasta.map { meta, _r, _a, _i, _fasta ->
            def dict = file(meta.genome_data.dict, checkIfExists: true)
            tuple(meta, dict)
        }

        ch_filtered_uBam
            .join(ch_cons_mapped_bam,    by: 0)
            .join(ch_cons_fasta_by_meta, by: 0)
            .join(ch_cons_dict_by_meta,  by: 0)
            .map { meta, ubam, mapped_bam, fasta, dict -> tuple(meta, ubam, mapped_bam, fasta, dict) }
            .set { ch_cons_zipperbam }

        FGBIO_ZIPPERBAMS_CONS(ch_cons_zipperbam)
        ch_versions = ch_versions.mix(FGBIO_ZIPPERBAMS_CONS.out.versions)

        FGBIO_SORTBAM(FGBIO_ZIPPERBAMS_CONS.out.bam)
        ch_versions = ch_versions.mix(FGBIO_SORTBAM.out.versions)

        def ch_consensus_filtered_bam = FGBIO_SORTBAM.out.bam

    // Consensus_filtered_bam into CRAM (integration into pipeline)

        SAMTOOLS_INDEX_CONSENSUS(ch_consensus_filtered_bam)
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX_CONSENSUS.out.versions)

        def ch_sam_convert_bai_fasta_fai = SAMTOOLS_INDEX_CONSENSUS.out.bai.map {meta, bai ->
            def gd    = (meta.genome_data instanceof Map) ? meta.genome_data : [:]
            def fasta = file(gd.fasta, checkIfExists: true)
            def fai   = file(gd.fai,   checkIfExists: true)
            tuple(meta, bai, fasta, fai)
        }

        def ch_consensus_bam_convert = ch_consensus_filtered_bam
            .join(ch_sam_convert_bai_fasta_fai, by: 0)

        SAMTOOLS_CONVERT(ch_consensus_bam_convert)
        ch_versions = ch_versions.mix(SAMTOOLS_CONVERT.out.versions)

        SAMTOOLS_CONVERT.out.cram
            .join(SAMTOOLS_CONVERT.out.crai, by: 0)
            .map { meta, cram, crai -> tuple(meta, cram, crai) }
            .set { ch_consensus_cram_crai }

    emit:
        ubam              = ch_ubam
        consensus_bam     = ch_consensus_filtered_bam
        grouped_bam       = ch_grouped_bam
        filtered_ubam     = ch_filtered_uBam
        consensus_cram_crai    = ch_consensus_cram_crai
        versions          = ch_versions
}
