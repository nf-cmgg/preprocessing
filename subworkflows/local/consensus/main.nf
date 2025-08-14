#!/usr/bin/env nextflow

include { FGBIO_COPYUMIFROMREADNAME         } from "${projectDir}/modules/nf-core/fgbio/copyumifromreadname/main.nf"
include { FGBIO_CALLDUPLEXCONSENSUSREADS    } from "${projectDir}/modules/nf-core/fgbio/callduplexconsensusreads/main.nf"
include { FGBIO_CALLMOLECULARCONSENSUSREADS } from "${projectDir}/modules/nf-core/fgbio/callmolecularconsensusreads/main.nf"
include { FGBIO_COLLECTDUPLEXSEQMETRICS     } from "${projectDir}/modules/nf-core/fgbio/collectduplexseqmetrics/main.nf"
include { FGBIO_FASTQTOBAM                  } from "${projectDir}/modules/nf-core/fgbio/fastqtobam/main.nf"
include { FGBIO_FILTERCONSENSUSREADS        } from "${projectDir}/modules/nf-core/fgbio/filterconsensusreads/main.nf"
include { FGBIO_GROUPREADSBYUMI             } from "${projectDir}/modules/nf-core/fgbio/groupreadsbyumi/main.nf"
include { FGBIO_SORTBAM                     } from "${projectDir}/modules/nf-core/fgbio/sortbam/main.nf"
include { FGBIO_ZIPPERBAMS                  } from "${projectDir}/modules/nf-core/fgbio/zipperbams/main.nf"

workflow CONSENSUS {
    take:
        ch_fastq                // channel: tuple(meta, fastq1, fastq2) for SE/PE/duplex samples
        ch_reference            // channel: reference genome fasta file
        ch_umi_in_readname      // boolean

    main:
        ch_versions       = Channel.empty()
        ch_duplex_metrics = Channel.empty()

        ch_fastq
          .combine(ch_umi_in_readname)
          .branch(
            rn:  { it[3] == true },
            seq: { it[3] == false }
          )
          .set { ch_fastq_branch }

        RN_FQ  = ch_fastq_branch.rn .map { meta, r1, r2, f -> [meta, r1, r2] }
        SEQ_FQ = ch_fastq_branch.seq.map { meta, r1, r2, f -> [meta, r1, r2] }

        FGBIO_FASTQTOBAM(RN_FQ, ch_reference)
        ch_versions = ch_versions.mix(FGBIO_FASTQTOBAM.out.versions.first())

        FGBIO_COPYUMIFROMREADNAME(FGBIO_FASTQTOBAM.out.bam)
        ch_versions = ch_versions.mix(FGBIO_COPYUMIFROMREADNAME.out.versions.first())

        FGBIO_SORTBAM(FGBIO_COPYUMIFROMREADNAME.out.bam)
        ch_versions = ch_versions.mix(FGBIO_SORTBAM.out.versions.first())

        FGBIO_GROUPREADSBYUMI(FGBIO_SORTBAM.out.bam)
        ch_versions = ch_versions.mix(FGBIO_GROUPREADSBYUMI.out.versions.first())

        FGBIO_CALLMOLECULARCONSENSUSREADS(FGBIO_GROUPREADSBYUMI.out.bam)
        ch_versions = ch_versions.mix(FGBIO_CALLMOLECULARCONSENSUSREADS.out.versions.first())

        if (params.enable_duplex) {
            FGBIO_ZIPPERBAMS(FGBIO_CALLMOLECULARCONSENSUSREADS.out.bam)
            ch_versions = ch_versions.mix(FGBIO_ZIPPERBAMS.out.versions.first())
            FGBIO_CALLDUPLEXCONSENSUSREADS(FGBIO_ZIPPERBAMS.out.bam)
            ch_versions = ch_versions.mix(FGBIO_CALLDUPLEXCONSENSUSREADS.out.versions.first())
            FGBIO_FILTERCONSENSUSREADS(FGBIO_CALLDUPLEXCONSENSUSREADS.out.bam)
            ch_versions = ch_versions.mix(FGBIO_FILTERCONSENSUSREADS.out.versions.first())
            FGBIO_COLLECTDUPLEXSEQMETRICS(FGBIO_FILTERCONSENSUSREADS.out.bam)
            ch_versions = ch_versions.mix(FGBIO_COLLECTDUPLEXSEQMETRICS.out.versions.first())
            ch_duplex_metrics = ch_duplex_metrics.mix(FGBIO_COLLECTDUPLEXSEQMETRICS.out.metrics)
        } else {
            FGBIO_FILTERCONSENSUSREADS(FGBIO_CALLMOLECULARCONSENSUSREADS.out.bam)
            ch_versions = ch_versions.mix(FGBIO_FILTERCONSENSUSREADS.out.versions.first())
        }

        ch_consensus_rn  = FGBIO_FILTERCONSENSUSREADS.out.bam
        ch_consensus_seq = Channel.empty()   // pending: ExtractUmisFromBam branch

        ch_consensus = ch_consensus_rn.mix(ch_consensus_seq)

    emit:
        consensus_bam  = ch_consensus
        duplex_metrics = ch_duplex_metrics
        versions       = ch_versions
}