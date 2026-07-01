#!/usr/bin/env nextflow

// MODULES
include { FGUMI_EXTRACT          } from "../../../modules/nf-core/fgumi/extract/main.nf"
include { FGUMI_FILTER           } from "../../../modules/nf-core/fgumi/filter/main.nf"
include { FGUMI_GROUP            } from "../../../modules/nf-core/fgumi/group/main.nf"
include { FGUMI_SIMPLEX          } from "../../../modules/nf-core/fgumi/simplex/main.nf"
include { CRAM_SNAPZIPPER_FGUMI  } from "../cram_snapzipper_fgumi/main.nf"

// FUNCTIONS
include { getGenomeAttribute      } from '../../local/utils_nfcore_preprocessing_pipeline'

workflow CRAM_UMICONSENSUS_FGUMI {
    take:
    ch_meta_reads_aligner_index_fasta // channel: [mandatory] [meta, reads, aligner, index, fasta]

    main:
    // Step numbers follow the fgumi basic workflow terminology (this path executes steps 1, 3, 4, 5, and 7).
    // Step 1: build an unmapped BAM with UMI tags from input FASTQ.
    FGUMI_EXTRACT(
        ch_meta_reads_aligner_index_fasta
            .map { meta, reads, _aligner, _index, _fasta, _fai -> [meta, reads, (meta.readgroup?.LB ?: meta.library ?: meta.id)] }
    )

    // Step 3: align with SNAP, zipper tags back, then template-coordinate sort.
    CRAM_SNAPZIPPER_FGUMI(
        FGUMI_EXTRACT.out.bam
            .join(
                ch_meta_reads_aligner_index_fasta.map { meta, _reads, _aligner, _index, fasta, fai ->
                    [meta, getGenomeAttribute(meta.genome_data, 'snap'), fasta, getGenomeAttribute(meta.genome_data, 'dict'), fai]
                },
            )
    )

    FGUMI_GROUP(
        CRAM_SNAPZIPPER_FGUMI.out.bam,
        (params.fgumi_group_strategy ?: 'adjacency')
    )

    FGUMI_SIMPLEX(
        FGUMI_GROUP.out.bam,
        (params.fgumi_simplex_min_reads ?: 1),
        false
    )

    // Step 7: filter consensus reads, then coordinate-sort/index for downstream CRAM conversion.
    FGUMI_FILTER(
        FGUMI_SIMPLEX.out.bam
            .join(ch_meta_reads_aligner_index_fasta)
            .map { meta, simplex_bams, _reads, _aligner, _index, fasta, _fai -> [meta, simplex_bams, fasta] },
        '1,1,1',
        false
    )

    emit:
    cram                  = FGUMI_FILTER.out.bam
    // Compatibility output kept for downstream interfaces; currently not produced by this branch.
    zipper_diagnostics    = channel.empty()
    grouping_metrics      = FGUMI_GROUP.out.metrics
    family_size_histogram = FGUMI_GROUP.out.histogram
    consensus_metrics     = FGUMI_SIMPLEX.out.stats
    filtering_metrics     = FGUMI_FILTER.out.stats
}
