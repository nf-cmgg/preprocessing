#!/usr/bin/env nextflow

// MODULES
include { FGUMI_EXTRACT          } from "../../../modules/nf-core/fgumi/extract/main.nf"
include { FGUMI_FILTER           } from "../../../modules/nf-core/fgumi/filter/main.nf"
include { FGUMI_GROUP            } from "../../../modules/nf-core/fgumi/group/main.nf"
include { FGUMI_SIMPLEX          } from "../../../modules/nf-core/fgumi/simplex/main.nf"
include { CRAM_SNAPZIPPER_FGUMI as RAW_CRAM_SNAPZIPPER_FGUMI  } from "../cram_snapzipper_fgumi/main.nf"
include { CRAM_SNAPZIPPER_FGUMI as UMI_CRAM_SNAPZIPPER_FGUMI  } from "../cram_snapzipper_fgumi/main.nf"

// FUNCTIONS
include { getGenomeAttribute      } from '../../local/utils_nfcore_preprocessing_pipeline'

workflow FASTQ_UMICONSENSUS_FGUMI {
    take:
    ch_meta_fastqs // channel: [mandatory] [meta, fastqs]

    main:
    // Step numbers follow the fgumi basic workflow terminology (this path executes steps 1, 3, 4, 5, and 7).
    // Step 1: build an unmapped BAM with UMI tags from input FASTQ.
    FGUMI_EXTRACT(
        ch_meta_fastqs
            .map { meta, fastqs -> [meta, fastqs, (meta.readgroup?.LB ?: meta.library ?: meta.id)] }
    )

    // Step 3: align with SNAP, zipper tags back, then template-coordinate sort.
    RAW_CRAM_SNAPZIPPER_FGUMI(
        FGUMI_EXTRACT.out.bam
            .join(ch_meta_fastqs)
            .map { meta, ubams, _fastqs ->
                [
                    meta,
                    ubams,
                    getGenomeAttribute(meta.genome_data, 'snap'),
                    getGenomeAttribute(meta.genome_data, 'fasta'),
                    getGenomeAttribute(meta.genome_data, 'dict'),
                    getGenomeAttribute(meta.genome_data, 'fai')
                ]
            },
    )

    FGUMI_GROUP(
        RAW_CRAM_SNAPZIPPER_FGUMI.out.bam,
        'adjacency'
    )

    FGUMI_SIMPLEX(
        FGUMI_GROUP.out.bam.map { meta, bams -> [ meta, bams, meta.fgumi_simplex_min_reads ] },
        false
    )

    // Step 7: filter consensus reads, then coordinate-sort/index for downstream CRAM conversion.
    FGUMI_FILTER(
        FGUMI_SIMPLEX.out.bam
            .join(ch_meta_fastqs)
            .map { meta, simplex_bams, _fastqs -> 
                [meta, simplex_bams, getGenomeAttribute(meta.genome_data, 'fasta')]
            },
        '1,1,1',
        false
    )

    UMI_CRAM_SNAPZIPPER_FGUMI(
        FGUMI_FILTER.out.bam
            .join(ch_meta_fastqs)
            .map { meta, filtered_bams, _fastqs -> 
                [
                    meta,
                    filtered_bams,
                    getGenomeAttribute(meta.genome_data, 'snap'),
                    getGenomeAttribute(meta.genome_data, 'fasta'),
                    getGenomeAttribute(meta.genome_data, 'dict'),
                    getGenomeAttribute(meta.genome_data, 'fai')
                ]
            }
    )

    emit:
    cram                  = UMI_CRAM_SNAPZIPPER_FGUMI.out.bam
    grouping_metrics      = FGUMI_GROUP.out.metrics
    family_size_histogram = FGUMI_GROUP.out.histogram
    consensus_metrics     = FGUMI_SIMPLEX.out.stats
    filtering_metrics     = FGUMI_FILTER.out.stats
}