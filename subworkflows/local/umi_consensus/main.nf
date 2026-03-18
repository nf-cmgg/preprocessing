include { FASTQ_ALIGN_DNA } from '../../nf-core/fastq_align_dna/main'
include { FASTQ_ALIGN_DNA as FASTQ_ALIGN_DNA_CONSENSUS } from '../../nf-core/fastq_align_dna/main'

include { SAMTOOLS_COLLATE as UMI_SAMTOOLS_COLLATE } from '../../../modules/nf-core/samtools/collate/main.nf'
include { SAMTOOLS_FIXMATE as UMI_SAMTOOLS_FIXMATE } from '../../../modules/nf-core/samtools/fixmate/main.nf'
include { SAMTOOLS_SORT as UMI_SAMTOOLS_SORT_TEMPLATE } from '../../../modules/nf-core/samtools/sort/main.nf'
include { SAMTOOLS_FASTQ as UMI_SAMTOOLS_FASTQ } from '../../../modules/nf-core/samtools/fastq/main.nf'
include { SAMTOOLS_SORT as UMI_SAMTOOLS_SORT_FINAL } from '../../../modules/nf-core/samtools/sort/main.nf'

include { FGBIO_COPYUMIFROMREADNAME as UMI_FGBIO_COPYUMIFROMREADNAME } from '../../../modules/nf-core/fgbio/copyumifromreadname/main.nf'
include { FGBIO_GROUPREADSBYUMI as UMI_FGBIO_GROUPREADSBYUMI } from '../../../modules/nf-core/fgbio/groupreadsbyumi/main.nf'
include { FGBIO_CALLMOLECULARCONSENSUSREADS as UMI_FGBIO_CALLMOLECULARCONSENSUSREADS } from '../../../modules/nf-core/fgbio/callmolecularconsensusreads/main.nf'
include { FGBIO_FILTERCONSENSUSREADS as UMI_FGBIO_FILTERCONSENSUSREADS } from '../../../modules/nf-core/fgbio/filterconsensusreads/main.nf'
include { FGBIO_ZIPPERBAMS as UMI_FGBIO_ZIPPERBAMS } from '../../../modules/nf-core/fgbio/zipperbams/main.nf'

include { UMI_LOCAL_SAMTOOLS_VIEW } from '../../../modules/local/umi_consensus/main.nf'

// UMI consensus workflow for DNA samples.
// Input channel shape:
//   [meta, reads, aligner, index, fasta]
// Output channels:
//   bam_bai      -> [meta, bam, bai]
//   family_sizes -> [meta, histogram]
workflow UMI_CONSENSUS_KAPA {
    take:
    ch_meta_reads_aligner_index_fasta // [meta, reads, aligner, index, fasta]

    main:
    // 1) Initial mapping of raw reads with the configured aligner/index.
    FASTQ_ALIGN_DNA(ch_meta_reads_aligner_index_fasta, false)

    // 2) Build reference helper channels reused by downstream modules.
    ch_meta_fasta_fai = ch_meta_reads_aligner_index_fasta
        .map { meta, _reads, _aligner, _index, fasta ->
            def fai = meta.genome_data?.fai ?: '/dev/null'
            [meta, fasta, file(fai, checkIfExists: true)]
        }

    ch_meta_fasta = ch_meta_reads_aligner_index_fasta
        .map { meta, _reads, _aligner, _index, fasta -> [meta, fasta] }

    ch_meta_dict = ch_meta_reads_aligner_index_fasta
        .map { meta, _reads, _aligner, _index, _fasta ->
            def dict = meta.genome_data?.dict ?: '/dev/null'
            [meta, file(dict, checkIfExists: true)]
        }

    // 3) Prepare read-pair metadata and UMI tags before consensus calling.
    UMI_LOCAL_SAMTOOLS_VIEW(FASTQ_ALIGN_DNA.out.bam)

    UMI_SAMTOOLS_COLLATE(
        UMI_LOCAL_SAMTOOLS_VIEW.out.bam,
        ch_meta_fasta_fai
    )

    UMI_SAMTOOLS_FIXMATE(UMI_SAMTOOLS_COLLATE.out.bam)

    UMI_SAMTOOLS_SORT_TEMPLATE(
        UMI_SAMTOOLS_FIXMATE.out.bam
            .join(ch_meta_fasta, by: 0)
            .map { meta, bam, fasta -> [meta, bam, fasta] },
        'bai'
    )

    ch_template_bam_bai = UMI_SAMTOOLS_SORT_TEMPLATE.out.bam
        .join(UMI_SAMTOOLS_SORT_TEMPLATE.out.bai, by: 0)
        .map { meta, bam, bai -> [meta, bam, bai] }

    // 4) Copy UMI from read names to RX tag, then group by UMI families.
    UMI_FGBIO_COPYUMIFROMREADNAME(ch_template_bam_bai)

    UMI_FGBIO_GROUPREADSBYUMI(
        UMI_FGBIO_COPYUMIFROMREADNAME.out.bam,
        channel.value('Adjacency')
    )

    // 5) Call and filter molecular consensus reads.
    UMI_FGBIO_CALLMOLECULARCONSENSUSREADS(
        UMI_FGBIO_GROUPREADSBYUMI.out.bam,
        UMI_FGBIO_GROUPREADSBYUMI.out.bam.map { meta, _bam -> meta.umi_min_reads ?: 2 },
        channel.value(20)
    )

    UMI_FGBIO_FILTERCONSENSUSREADS(
        UMI_FGBIO_CALLMOLECULARCONSENSUSREADS.out.bam,
        ch_meta_fasta,
        UMI_FGBIO_CALLMOLECULARCONSENSUSREADS.out.bam.map { meta, _bam -> meta.umi_min_reads ?: 2 },
        channel.value(45),
        channel.value(0.2)
    )

    UMI_SAMTOOLS_FASTQ(
        UMI_FGBIO_FILTERCONSENSUSREADS.out.bam,
        true
    )

    // 6) Re-map consensus reads with the same aligner/index used upstream.
    FASTQ_ALIGN_DNA_CONSENSUS(
        UMI_SAMTOOLS_FASTQ.out.interleaved
            .join(ch_meta_reads_aligner_index_fasta.map { meta, _reads, aligner, index, fasta -> [meta, aligner, index, fasta] }, by: 0)
            .map { meta, interleaved_fastq, aligner, index, fasta -> [meta, [interleaved_fastq], aligner, index, fasta] }
    , false)

    FASTQ_ALIGN_DNA_CONSENSUS.out.bam
        .join(UMI_FGBIO_FILTERCONSENSUSREADS.out.bam, by: 0)
        .join(ch_meta_fasta, by: 0)
        .join(ch_meta_dict, by: 0)
        .map { meta, remap_bam, unmapped_bam, fasta, dict -> [meta, unmapped_bam, remap_bam, fasta, dict] }
        .multiMap { meta, unmapped_bam, remap_bam, fasta, dict ->
            unmapped: [meta, unmapped_bam]
            mapped: [meta, remap_bam]
            fasta: [meta, fasta]
            dict: [meta, dict]
        }
        .set { ch_zipper_inputs }

    // 7) Transfer unmapped metadata back to mapped consensus alignments.
    UMI_FGBIO_ZIPPERBAMS(
        ch_zipper_inputs.unmapped,
        ch_zipper_inputs.mapped,
        ch_zipper_inputs.fasta,
        ch_zipper_inputs.dict
    )

    // 8) Final coordinate sort + index for downstream CRAM conversion.
    UMI_SAMTOOLS_SORT_FINAL(
        UMI_FGBIO_ZIPPERBAMS.out.bam
            .join(ch_meta_fasta, by: 0)
            .map { meta, bam, fasta -> [meta, bam, fasta] },
        'bai'
    )

    emit:
    bam_bai = UMI_SAMTOOLS_SORT_FINAL.out.bam
        .join(UMI_SAMTOOLS_SORT_FINAL.out.bai, by: 0)
    family_sizes = UMI_FGBIO_GROUPREADSBYUMI.out.histogram
}
