#!/usr/bin/env nextflow

// MODULES
include { FGUMI_EXTRACT          } from "../../../modules/nf-core/fgumi/extract/main.nf"
include { FGUMI_FILTER           } from "../../../modules/nf-core/fgumi/filter/main.nf"
include { FGUMI_GROUP            } from "../../../modules/nf-core/fgumi/group/main.nf"
include { FGUMI_SIMPLEX          } from "../../../modules/nf-core/fgumi/simplex/main.nf"
include { FGUMI_SNAP_ZIPPER      } from "../../../modules/local/fgumi/snapzipper/main.nf"
include { SAMTOOLS_SORT          } from "../../../modules/nf-core/samtools/sort/main.nf"

// FUNCTIONS
include { getGenomeAttribute      } from '../../local/utils_nfcore_preprocessing_pipeline'

workflow UMI_CONSENSUS_FGUMI {
    take:
    ch_meta_reads_aligner_index_fasta // channel: [mandatory] [meta, reads, aligner, index, fasta]

    main:
    // Step 1: build an unmapped BAM with UMI tags from input FASTQ.
    FGUMI_EXTRACT(
        ch_meta_reads_aligner_index_fasta
            .map { meta, reads, _aligner, _index, _fasta -> [meta, reads, (meta.readgroup?.LB ?: meta.library ?: meta.id)] }
    )

    // Step 3: align with SNAP, zipper tags back, then template-coordinate sort.
    FGUMI_SNAP_ZIPPER(
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
        FGUMI_SNAP_ZIPPER.out.bam,
        (params.fgumi_group_strategy ?: 'adjacency')
    )

    FGUMI_SIMPLEX(
        FGUMI_GROUP.out.bam,
        (params.fgumi_simplex_min_reads ?: 1),
        false
    )

    // Step 7: filter consensus reads, then coordinate-sort/index for downstream CRAM conversion.
    FGUMI_FILTER(
        FGUMI_SIMPLEX.out.bam,
        FGUMI_SIMPLEX.out.bam
            .join(
                ch_meta_reads_aligner_index_fasta.map { meta, _reads, _aligner, _index, fasta -> [meta, fasta] },
                by: 0,
            )
            .map { meta, _bam, fasta -> [meta, fasta] },
        "1,1,1",
        false
    )

    SAMTOOLS_SORT(
        FGUMI_FILTER.out.bam
            .join(
                ch_meta_reads_aligner_index_fasta.map { meta, _reads, _aligner, _index, fasta -> [meta, fasta] },
                by: 0,
            )
            .map { meta, bam, fasta -> [meta, bam, fasta] },
        "crai"
    )

    emit:
    cram_crai             = SAMTOOLS_SORT.out.cram.join(SAMTOOLS_SORT.out.crai, failOnMismatch: true, failOnDuplicate: true)
    grouping_metrics      = FGUMI_GROUP.out.metrics
    family_size_histogram = FGUMI_GROUP.out.histogram
    consensus_metrics     = FGUMI_SIMPLEX.out.stats
    filtering_metrics     = FGUMI_FILTER.out.stats
    filtered_consensus_cram = SAMTOOLS_SORT.out.cram
}
