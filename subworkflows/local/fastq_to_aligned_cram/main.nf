#!/usr/bin/env nextflow

//
// Take fastq; align, postprocess and compress
//

// MODULES
include { BIOBAMBAM_BAMSORMADUP } from "../../../modules/nf-core/biobambam/bamsormadup/main.nf"
include { SAMTOOLS_CONVERT      } from "../../../modules/nf-core/samtools/convert/main"
include { SAMTOOLS_SORMADUP     } from "../../../modules/nf-core/samtools/sormadup/main.nf"
include { SAMTOOLS_SORT         } from "../../../modules/nf-core/samtools/sort/main"

// SUBWORKFLOWS
include { FASTQ_ALIGN_DNA       } from '../../nf-core/fastq_align_dna/main'
include { FASTQ_ALIGN_RNA       } from '../../local/fastq_align_rna/main'
include { CRAM_UMICONSENSUS_FGUMI   } from '../../local/cram_umiconsensus_fgumi/main.nf'

// FUNCTIONS
include { getGenomeAttribute    } from '../../local/utils_nfcore_preprocessing_pipeline'

workflow FASTQ_TO_CRAM {
    take:
    ch_meta_reads_aligner_index_fasta_fai_gtf // channel: [mandatory] [meta, [fastq, ...], aligner [bowtie2, bwamem, bwamem2, dragmap, snap, star], aligner_index, fasta, fai, gtf]

    main:
    ch_sormadup_metrics = channel.empty()

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // STEP: ALIGNMENT
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_meta_reads_aligner_index_fasta_fai_gtf.dump(tag: "FASTQ_TO_CRAM: reads to align", pretty: true)
    ch_meta_reads_aligner_index_fasta_fai_gtf
        .branch { meta, reads, aligner, index, fasta, fai, gtf ->
            rna: meta.sample_type == "RNA"
            return [meta, reads, "star", getGenomeAttribute(meta.genome_data, 'star'), gtf]
            dna: true
            // catch all non-RNA samples as DNA, as some may be missing sample_type or have other sample types (e.g. tissue, cell line, etc.) that should be aligned with the DNA aligner
            //dna: meta.sample_type == "DNA" || meta.sample_type == "Tissue"
            return [meta, reads, aligner, index, fasta, fai]
        }
        .set { ch_meta_reads_aligner_index_fasta_datatype }

    ch_meta_reads_aligner_index_fasta_datatype.dna
        .branch { meta, _reads, _aligner, _index, _fasta, _fai ->
            // fgumi consensus is opt-in via fgumi_aware to avoid changing samtools umi_aware semantics.
            umi: meta.fgumi_aware
            non_umi: true
        }
        .set { ch_dna_to_align }

    // Align non-UMI DNA fastq files per sample
    FASTQ_ALIGN_DNA(
        ch_dna_to_align.non_umi,
        false,
    )

    // UMI-aware fgumi branch (steps 1, 3, 4, 5, 6, 7 in fgumi Basic Workflow)
    CRAM_UMICONSENSUS_FGUMI(
        ch_dna_to_align.umi
    )

    FASTQ_ALIGN_RNA(
        ch_meta_reads_aligner_index_fasta_datatype.rna
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // STEP: MARK DUPLICATES
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    FASTQ_ALIGN_DNA.out.bam
        .mix(FASTQ_ALIGN_RNA.out.bam)
        .map { meta, files ->
            def gk = (meta.chunks as Integer ?: 1)
            return [
                groupKey(
                    meta - meta.subMap('readgroup', 'chunks') + [id: meta.id ==~ /^\d{4}\..*$/ ? meta.id[5..-1] : meta.id],
                    gk,
                ),
                files,
            ]
        }
        .groupTuple()
        .map { meta, files ->
            def gk = (meta.count as Integer ?: 1)
            return [
                groupKey(
                    meta - meta.subMap('count') + [id: meta.samplename ?: meta.id],
                    gk,
                ),
                files,
            ]
        }
        .groupTuple()
        .map { meta, files ->
            return [meta, files.flatten(), getGenomeAttribute(meta.genome_data, 'fasta'), getGenomeAttribute(meta.genome_data, 'fai')]
        }
        .dump(tag: "FASTQ_TO_CRAM: aligned bam per sample", pretty: true)
        .branch { meta, files, fasta, fai ->
            bamsormadup: meta.markdup == "bamsormadup"
            return [meta, files, fasta, fai]
            samtools: meta.markdup == "samtools"
            return [meta, files, fasta, fai]
            sort: meta.markdup == "false" || meta.markdup == false
            return [meta, files, fasta, fai]
            unknown: true
            error("markdup option ${meta.markdup} not supported")
        }
        .set { ch_bam_fasta_fai }

    ch_markdup_index = channel.empty()

    // UMI branch outputs are mixed into the common markdup/metrics streams.
    ch_markdup_index = ch_markdup_index.mix(
        CRAM_UMICONSENSUS_FGUMI.out.cram_crai
    )
    ch_sormadup_metrics = ch_sormadup_metrics.mix(CRAM_UMICONSENSUS_FGUMI.out.grouping_metrics)
    ch_sormadup_metrics = ch_sormadup_metrics.mix(CRAM_UMICONSENSUS_FGUMI.out.family_size_histogram)
    ch_sormadup_metrics = ch_sormadup_metrics.mix(CRAM_UMICONSENSUS_FGUMI.out.consensus_metrics)
    ch_sormadup_metrics = ch_sormadup_metrics.mix(CRAM_UMICONSENSUS_FGUMI.out.filtering_metrics)
    ch_family_size_histogram = CRAM_UMICONSENSUS_FGUMI.out.family_size_histogram
    ch_filtered_consensus_cram = CRAM_UMICONSENSUS_FGUMI.out.filtered_consensus_cram
    ch_zipper_diagnostics = CRAM_UMICONSENSUS_FGUMI.out.zipper_diagnostics

    // BIOBAMBAM_BAMSORMADUP([meta, [bam, bam]], fasta, fai)
    BIOBAMBAM_BAMSORMADUP(ch_bam_fasta_fai.bamsormadup)
    ch_markdup_index = ch_markdup_index.mix(BIOBAMBAM_BAMSORMADUP.out.bam.join(BIOBAMBAM_BAMSORMADUP.out.bam_index, failOnMismatch: true, failOnDuplicate: true))
    ch_sormadup_metrics = ch_sormadup_metrics.mix(BIOBAMBAM_BAMSORMADUP.out.metrics)

    // SAMTOOLS_SORMADUP([meta, [bam, bam]], fasta, fai)
    SAMTOOLS_SORMADUP(ch_bam_fasta_fai.samtools)
    ch_markdup_index = ch_markdup_index.mix(SAMTOOLS_SORMADUP.out.cram.join(SAMTOOLS_SORMADUP.out.crai, failOnMismatch: true, failOnDuplicate: true))
    ch_sormadup_metrics = ch_sormadup_metrics.mix(SAMTOOLS_SORMADUP.out.metrics)

    // Merge bam files and compress
    // SAMTOOLS_SORT([meta, [bam, bam], fasta],index_format)
    SAMTOOLS_SORT(ch_bam_fasta_fai.sort, "crai")
    ch_markdup_index = ch_markdup_index.mix(SAMTOOLS_SORT.out.cram.join(SAMTOOLS_SORT.out.index, failOnMismatch: true, failOnDuplicate: true))

    ch_markdup_index.dump(tag: "FASTQ_TO_CRAM: postprocessed bam", pretty: true)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // COMPRESSION
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_markdup_index
        .branch { meta, reads, index ->
            bam: reads.getExtension() == "bam"
            return [meta, reads, index]
            cram: reads.getExtension() == "cram"
            return [meta, reads, index]
        }
        .set { ch_markdup_index }

    ch_markdup_index.bam
        .map { meta, bam, bai ->
            bam_bai: [meta, bam, bai, getGenomeAttribute(meta.genome_data, 'fasta'), getGenomeAttribute(meta.genome_data, 'fai')]
        }
        .set { ch_bam_bai_fasta_fai }

    SAMTOOLS_CONVERT(ch_bam_bai_fasta_fai)

    ch_markdup_index.cram
        .mix(
            SAMTOOLS_CONVERT.out.cram.join(SAMTOOLS_CONVERT.out.crai, failOnMismatch: true, failOnDuplicate: true)
        )
        .set { ch_cram_crai }
    ch_cram_crai.dump(tag: "FASTQ_TO_CRAM: cram and crai", pretty: true)

    emit:
    cram_crai            = ch_cram_crai
    // UMI-specific output channels for downstream reporting and publishing.
    filtered_consensus_cram = ch_filtered_consensus_cram
    zipper_diagnostics   = ch_zipper_diagnostics
    rna_splice_junctions = FASTQ_ALIGN_RNA.out.splice_junctions
    rna_junctions        = FASTQ_ALIGN_RNA.out.junctions
    sormadup_metrics     = ch_sormadup_metrics
    family_size_histogram = ch_family_size_histogram
    align_reports        = FASTQ_ALIGN_DNA.out.reports
}
