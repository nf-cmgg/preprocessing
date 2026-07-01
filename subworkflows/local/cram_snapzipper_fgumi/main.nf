#!/usr/bin/env nextflow

// MODULES
include { FGUMI_SNAP_ALIGN                   } from "../../../modules/local/fgumi/snapalign/main.nf"
include { FGUMI_ZIPPER                        } from "../../../modules/local/fgumi/zipper/main.nf"
include { FGUMI_SORT as FGUMI_TEMPLATE_SORT  } from "../../../modules/nf-core/fgumi/sort/main.nf"
include { SAMTOOLS_SORT as SAMTOOLS_QNAME_SORT_UNMAPPED } from "../../../modules/nf-core/samtools/sort/main.nf"
include { SAMTOOLS_SORT as SAMTOOLS_QNAME_SORT_MAPPED   } from "../../../modules/nf-core/samtools/sort/main.nf"

workflow CRAM_SNAPZIPPER_FGUMI {
    take:
    ch_meta_unmapped_index_fasta_dict_fai

    main:
    FGUMI_SNAP_ALIGN(ch_meta_unmapped_index_fasta_dict_fai.map { meta, unmapped_bam, index, fasta, dict, _fai -> [meta, unmapped_bam, index, fasta, dict] })

    // Queryname sort the unmapped BAM in parallel with mapped BAM sort.
    SAMTOOLS_QNAME_SORT_UNMAPPED(
        ch_meta_unmapped_index_fasta_dict_fai.map { meta, unmapped_bam, _index, fasta, _dict, fai -> [meta, unmapped_bam, fasta, fai] },
        ''
    )

    // Sort mapped alignments by queryname and emit SAM for zipper stdin.
    SAMTOOLS_QNAME_SORT_MAPPED(
        FGUMI_SNAP_ALIGN.out.mapped_bam
            .join(
                ch_meta_unmapped_index_fasta_dict_fai.map { meta, _unmapped_bam, _index, fasta, _dict, fai -> [meta, fasta, fai] },
            ),
        ''
    )

    FGUMI_ZIPPER(
        SAMTOOLS_QNAME_SORT_MAPPED.out.bam
            .join(SAMTOOLS_QNAME_SORT_UNMAPPED.out.bam)
            .join(
                ch_meta_unmapped_index_fasta_dict_fai.map { meta, _unmapped_bam, _index, fasta, dict, _fai -> [meta, fasta, dict] },
            )
    )

    FGUMI_TEMPLATE_SORT(FGUMI_ZIPPER.out.bam)

    emit:
    bam            = FGUMI_TEMPLATE_SORT.out.bam
    versions_fgumi = FGUMI_SNAP_ALIGN.out.versions_fgumi.mix(FGUMI_ZIPPER.out.versions_fgumi).mix(FGUMI_TEMPLATE_SORT.out.versions_fgumi)
}
