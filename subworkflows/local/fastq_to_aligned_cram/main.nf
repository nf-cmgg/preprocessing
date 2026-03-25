#!/usr/bin/env nextflow

//
// Take fastq; align, postprocess and compress
//

// MODULES
include { BIOBAMBAM_BAMSORMADUP } from "../../../modules/nf-core/biobambam/bamsormadup/main.nf"
include { BWA_MEM as FASTQ_ALIGN_DNA_BWAMEM } from '../../../modules/nf-core/bwa/mem/main.nf'
include { SNAPALIGNER_ALIGN as SNAP_ALIGN } from '../../../modules/local/snapaligner/align/main.nf'
include { SAMTOOLS_CONVERT      } from "../../../modules/nf-core/samtools/convert/main"
include { SAMTOOLS_SORMADUP     } from "../../../modules/nf-core/samtools/sormadup/main.nf"
include { SAMTOOLS_SORT         } from "../../../modules/nf-core/samtools/sort/main"
include { SAMTOOLS_SORT as UMI_SAMTOOLS_SORT_PREP } from "../../../modules/nf-core/samtools/sort/main"
include { SAMTOOLS_SORT as UMI_SAMTOOLS_MERGE_CRAM } from "../../../modules/nf-core/samtools/sort/main"
include { UMI_CONSENSUS_KAPA    } from '../../local/umi_consensus/main'

// SUBWORKFLOWS
include { FASTQ_ALIGN_RNA       } from '../../local/fastq_align_rna/main'

// FUNCTIONS
include { getGenomeAttribute    } from '../../local/utils_nfcore_preprocessing_pipeline'

workflow FASTQ_ALIGN_DNA {
    take:
    ch_meta_reads_aligner_index_fasta // [meta, reads, aligner, index, fasta]
    sort

    main:
    // Validate aligner values and normalize into per-aligner tuple shapes.
    ch_meta_reads_aligner_index_fasta
        .map { meta, reads, aligner, index, fasta ->
            if (!(aligner in ['bwamem', 'snap'])) {
                error("FASTQ_ALIGN_DNA currently supports aligners 'bwamem' and 'snap', got: ${aligner}")
            }
            [meta, reads, aligner, index, fasta]
        }
        .branch { meta, reads, aligner, index, fasta ->
            bwamem: aligner == 'bwamem'
            return [meta, reads, index, fasta]
            snap: aligner == 'snap'
            return [meta, reads, index]
        }
        .set { ch_align }

    // Run aligners independently and merge back into one output contract.
    FASTQ_ALIGN_DNA_BWAMEM(ch_align.bwamem, sort)
    SNAP_ALIGN(ch_align.snap)

    emit:
    bam = FASTQ_ALIGN_DNA_BWAMEM.out.bam.mix(SNAP_ALIGN.out.bam)
    bam_index = FASTQ_ALIGN_DNA_BWAMEM.out.csi.mix(SNAP_ALIGN.out.bai)
    reports = channel.empty()
}

workflow FASTQ_TO_CRAM {
    take:
    ch_meta_reads_aligner_index_fasta_gtf // channel: [mandatory] [meta, [fastq, ...], aligner [bowtie2, bwamem, bwamem2, dragmap, snap, star], aligner_index, fasta, gtf]

    main:
    ch_sormadup_metrics = channel.empty()
    ch_umi_family_sizes = channel.empty()

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // STEP: ALIGNMENT
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Split by analysis path:
    // - umi: DNA samples with UMI consensus enabled
    // - rna: RNA samples (forced to STAR)
    // - dna: all remaining non-RNA samples
    ch_meta_reads_aligner_index_fasta_gtf.dump(tag: "FASTQ_TO_CRAM: reads to align", pretty: true)
    ch_meta_reads_aligner_index_fasta_gtf
        .branch { meta, reads, aligner, index, fasta, gtf ->
            umi: meta.umi_consensus && meta.sample_type != "RNA"
            return [meta, reads, aligner, index, fasta]
            rna: meta.sample_type == "RNA"
            return [meta, reads, "star", getGenomeAttribute(meta.genome_data, 'star'), gtf]
            dna: true
            // catch all non-RNA samples as DNA, as some may be missing sample_type or have other sample types (e.g. tissue, cell line, etc.) that should be aligned with the DNA aligner
            //dna: meta.sample_type == "DNA" || meta.sample_type == "Tissue"
            return [meta, reads, aligner, index, fasta]
        }
        .set { ch_meta_reads_aligner_index_fasta_datatype }

    // Align non-RNA samples with DNA aligners and RNA samples with STAR.
    ch_dna_umi = ch_meta_reads_aligner_index_fasta_datatype.dna.mix(ch_meta_reads_aligner_index_fasta_datatype.umi)

    FASTQ_ALIGN_DNA(
        ch_dna_umi,
        false,
    )
    FASTQ_ALIGN_RNA(
        ch_meta_reads_aligner_index_fasta_datatype.rna
    )

    FASTQ_ALIGN_DNA.out.bam
        .filter { meta, _bam -> meta.umi_consensus && meta.sample_type != "RNA" }
        .join(ch_meta_reads_aligner_index_fasta_datatype.umi.map { meta, _reads, _aligner, _index, fasta -> [meta, fasta] }, by: 0)
        .map { meta, bam, fasta -> [meta, bam, fasta] }
        .set { ch_umi_bam_fasta }

    // UMI path pre-sort: consensus input BAM must be coordinate-sorted/indexed.
    UMI_SAMTOOLS_SORT_PREP(
        ch_umi_bam_fasta,
        'bai'
    )

    // Build UMI consensus per chunk (preserving aligner + index used upstream).
    UMI_CONSENSUS_KAPA(
        UMI_SAMTOOLS_SORT_PREP.out.bam
            .join(UMI_SAMTOOLS_SORT_PREP.out.bai, by: 0)
            .join(ch_meta_reads_aligner_index_fasta_datatype.umi.map { meta, _reads, aligner, index, fasta -> [meta, aligner, index, fasta] }, by: 0)
            .map { meta, bam, bai, aligner, index, fasta -> [meta, bam, bai, aligner, index, fasta] }
    )

    // Chunk-level UMI BAM outputs are converted to CRAM later via SAMTOOLS_CONVERT.
    UMI_CONSENSUS_KAPA.out.bam_bai
        .map { meta, bam, bai ->
            [meta, bam, bai, getGenomeAttribute(meta.genome_data, 'fasta'), getGenomeAttribute(meta.genome_data, 'fai')]
        }
        .set { ch_umi_bam_bai_fasta_fai }

    // Merge chunk BAMs into one per-sample CRAM for UMI samples.
    // Group key is sample-level id (`samplename` fallback `id`).
    UMI_CONSENSUS_KAPA.out.bam_bai
        .map { meta, bam, _bai ->
            [meta.samplename ?: meta.id, meta, bam, getGenomeAttribute(meta.genome_data, 'fasta')]
        }
        .groupTuple(by: 0)
        .map { _sample_id, metas, bams, fastas ->
            def base_meta = metas[0]
            def merged_meta = base_meta - base_meta.subMap('readgroup', 'chunks') + [id: base_meta.samplename ?: base_meta.id]
            [merged_meta, bams.flatten(), fastas[0]]
        }
        .set { ch_umi_merge_input }

    UMI_SAMTOOLS_MERGE_CRAM(
        ch_umi_merge_input,
        'crai'
    )

    UMI_SAMTOOLS_MERGE_CRAM.out.cram
        .join(UMI_SAMTOOLS_MERGE_CRAM.out.crai, by: 0)
        .set { ch_umi_cram_crai_merged }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // STEP: MARK DUPLICATES
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Non-UMI path for duplicate marking / direct merge-sort handling.
    FASTQ_ALIGN_DNA.out.bam
        .filter { meta, _files -> !(meta.umi_consensus && meta.sample_type != "RNA") }
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
        .set { ch_bam_fasta }

    ch_markdup_index = channel.empty()

    // markdup=bamsormadup
    BIOBAMBAM_BAMSORMADUP(ch_bam_fasta.bamsormadup)
    ch_markdup_index = ch_markdup_index.mix(BIOBAMBAM_BAMSORMADUP.out.bam.join(BIOBAMBAM_BAMSORMADUP.out.bam_index, failOnMismatch: true, failOnDuplicate: true))
    ch_sormadup_metrics = ch_sormadup_metrics.mix(BIOBAMBAM_BAMSORMADUP.out.metrics)

    // markdup=samtools
    SAMTOOLS_SORMADUP(ch_bam_fasta.samtools)
    ch_markdup_index = ch_markdup_index.mix(SAMTOOLS_SORMADUP.out.cram.join(SAMTOOLS_SORMADUP.out.crai, failOnMismatch: true, failOnDuplicate: true))
    ch_sormadup_metrics = ch_sormadup_metrics.mix(SAMTOOLS_SORMADUP.out.metrics)

    // markdup=false: merge/sort only (no duplicate marking)
    SAMTOOLS_SORT(ch_bam_fasta.sort, "crai")
    ch_markdup_index = ch_markdup_index.mix(SAMTOOLS_SORT.out.cram.join(SAMTOOLS_SORT.out.crai, failOnMismatch: true, failOnDuplicate: true))

    ch_markdup_index.dump(tag: "FASTQ_TO_CRAM: postprocessed bam", pretty: true)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // COMPRESSION
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Normalize outputs into bam/cram branches before final CRAM aggregation.
    ch_markdup_index
        .branch { meta, reads, index ->
            bam: reads.getExtension() == "bam"
            return [meta, reads, index]
            cram: reads.getExtension() == "cram"
            return [meta, reads, index]
        }
        .set { ch_markdup_index }

    // Convert remaining BAM outputs to CRAM (including UMI chunk BAMs).
    ch_markdup_index.bam
        .map { meta, bam, bai ->
            bam_bai: [meta, bam, bai, getGenomeAttribute(meta.genome_data, 'fasta'), getGenomeAttribute(meta.genome_data, 'fai')]
        }
        .mix(ch_umi_bam_bai_fasta_fai)
        .set { ch_bam_bai_fasta_fai }

    SAMTOOLS_CONVERT(ch_bam_bai_fasta_fai)

    ch_converted_cram_crai = SAMTOOLS_CONVERT.out.cram
        .join(SAMTOOLS_CONVERT.out.crai, failOnMismatch: true, failOnDuplicate: true)

    // Keep UMI chunk CRAMs separate for explicit publication.
    ch_umi_cram_crai_chunks = ch_converted_cram_crai
        .filter { meta, _cram, _crai -> meta.umi_consensus && meta.sample_type != "RNA" }

    // Non-UMI CRAMs from converted BAM path.
    ch_non_umi_cram_crai = ch_converted_cram_crai
        .filter { meta, _cram, _crai -> !(meta.umi_consensus && meta.sample_type != "RNA") }

    // Final main CRAM channel contains:
    // - CRAMs emitted directly by markdup/sort paths
    // - non-UMI CRAMs converted from BAM
    // - merged per-sample UMI CRAMs
    ch_markdup_index.cram
        .mix(ch_non_umi_cram_crai)
        .mix(ch_umi_cram_crai_merged)
        .set { ch_cram_crai }
    ch_cram_crai.dump(tag: "FASTQ_TO_CRAM: cram and crai", pretty: true)

    ch_umi_family_sizes = ch_umi_family_sizes.mix(UMI_CONSENSUS_KAPA.out.family_sizes)

    emit:
    cram_crai            = ch_cram_crai
    umi_cram_crai_chunks = ch_umi_cram_crai_chunks
    umi_cram_crai_merged = ch_umi_cram_crai_merged
    rna_splice_junctions = FASTQ_ALIGN_RNA.out.splice_junctions
    rna_junctions        = FASTQ_ALIGN_RNA.out.junctions
    sormadup_metrics     = ch_sormadup_metrics
    umi_family_sizes     = ch_umi_family_sizes
    align_reports        = FASTQ_ALIGN_DNA.out.reports
}
