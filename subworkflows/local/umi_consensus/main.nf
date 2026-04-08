#!/usr/bin/env nextflow

// MODULES
include { FGUMI_DUPLEX_METRICS   } from "../../../modules/local/fgumi/duplexmetrics/main.nf"
include { FGUMI_EXTRACT          } from "../../../modules/local/fgumi/extract/main.nf"
include { FGUMI_FILTER           } from "../../../modules/local/fgumi/filter/main.nf"
include { FGUMI_GROUP            } from "../../../modules/local/fgumi/group/main.nf"
include { FGUMI_SIMPLEX          } from "../../../modules/local/fgumi/simplex/main.nf"
include { FGUMI_SNAP_ZIPPER_SORT } from "../../../modules/local/fgumi/snapzippersort/main.nf"
include { FGUMI_SORT             } from "../../../modules/local/fgumi/sort/main.nf"

// FUNCTIONS
include { getGenomeAttribute      } from '../../local/utils_nfcore_preprocessing_pipeline'

workflow UMI_CONSENSUS_FGUMI {
    take:
    ch_meta_reads_aligner_index_fasta // channel: [mandatory] [meta, reads, aligner, index, fasta]

    main:
    // Step 1: build an unmapped BAM with UMI tags from input FASTQ.
    FGUMI_EXTRACT(
        ch_meta_reads_aligner_index_fasta
            .map { meta, reads, _aligner, _index, _fasta -> [meta, reads] }
    )

    // Step 3: align with SNAP, zipper tags back, then template-coordinate sort.
    FGUMI_SNAP_ZIPPER_SORT(
        FGUMI_EXTRACT.out.bam
            .join(
                ch_meta_reads_aligner_index_fasta.map { meta, _reads, _aligner, _index, fasta ->
                    [meta, getGenomeAttribute(meta.genome_data, 'snap'), fasta, getGenomeAttribute(meta.genome_data, 'dict')]
                },
                by: 0,
            )
            .map { meta, unmapped_bam, snap_index, fasta, dict -> [meta, unmapped_bam, snap_index, fasta, dict] }
    )

    FGUMI_GROUP(
        FGUMI_SNAP_ZIPPER_SORT.out.bam
    )

    FGUMI_SIMPLEX(
        FGUMI_GROUP.out.bam
    )

    FGUMI_DUPLEX_METRICS(
        FGUMI_GROUP.out.bam
    )

    // Step 7: filter consensus reads, then coordinate-sort/index for downstream CRAM conversion.
    FGUMI_FILTER(
        FGUMI_SIMPLEX.out.bam
            .join(
                ch_meta_reads_aligner_index_fasta.map { meta, _reads, _aligner, _index, fasta -> [meta, fasta] },
                by: 0,
            )
            .map { meta, bam, fasta -> [meta, bam, fasta] }
    )

    FGUMI_SORT(
        FGUMI_FILTER.out.bam
    )

    emit:
    bam_bai               = FGUMI_SORT.out.bam.join(FGUMI_SORT.out.bai, failOnMismatch: true, failOnDuplicate: true)
    grouping_metrics      = FGUMI_GROUP.out.grouping_metrics
    family_size_histogram = FGUMI_GROUP.out.family_size_histogram
    consensus_metrics     = FGUMI_SIMPLEX.out.consensus_metrics
    filtering_metrics     = FGUMI_FILTER.out.filtering_metrics
    duplex_metrics        = FGUMI_DUPLEX_METRICS.out.duplex_metrics
    filtered_consensus_bam = FGUMI_SORT.out.bam
}
