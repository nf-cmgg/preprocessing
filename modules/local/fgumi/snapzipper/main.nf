include { FGUMI_SORT as FGUMI_TEMPLATE_SORT } from "../../../nf-core/fgumi/sort/main.nf"
include { FGUMI_ZIPPER                         } from "../zipper/main.nf"
include { SAMTOOLS_SORT as SAMTOOLS_QNAME_SORT_UNMAPPED } from "../../../nf-core/samtools/sort/main.nf"
include { SAMTOOLS_SORT as SAMTOOLS_QNAME_SORT_MAPPED   } from "../../../nf-core/samtools/sort/main.nf"

process FGUMI_SNAP_ZIPPER_RUN {
    tag "$meta.id"
    label 'process_high'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/763a833519c23555be888065f492215f57344155106972e272a0f8df78c57659/data'
        : 'community.wave.seqera.io/library/fgumi_samtools_snap-aligner:c985f9394623a414'}"

    input:
    tuple val(meta), path(unmapped_bam), path(index, stageAs: "index/*"), path(fasta), path(dict)

    output:
    tuple val(meta), path("${prefix}.snap.bam"), emit: mapped_bam
    tuple val(meta), path(unmapped_bam), emit: unmapped_bam
    tuple val("${task.process}"), val('fgumi'), eval("fgumi --version | sed 's/^fgumi //;q'"), topic: versions, emit: versions_fgumi

    when:
    task.ext.when == null || task.ext.when

    script:
    def snap_args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.fgumi"

    """
    INDEX_FILE=\$(find -L ./ -name "OverflowTable*" -print -quit)
    [ -z "\$INDEX_FILE" ] && echo "Snap index files not found" 1>&2 && exit 1
    INDEX=\$(dirname "\$INDEX_FILE")

    fgumi fastq --input ${unmapped_bam} \
        | snap-aligner paired \
            \$INDEX \
            -pairedInterleavedFastq - \
            -o -sam - \
            -t ${task.cpus} \
            ${snap_args} \
        > ${prefix}.snap.sam

    samtools view \
        -@ ${task.cpus} \
        -b \
        -o ${prefix}.snap.bam \
        ${prefix}.snap.sam
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.fgumi"
    """
    touch ${prefix}.snap.bam
    touch ${unmapped_bam}
    """
}

workflow FGUMI_SNAP_ZIPPER {
    take:
    ch_meta_unmapped_index_fasta_dict

    main:
    FGUMI_SNAP_ZIPPER_RUN(ch_meta_unmapped_index_fasta_dict)

    // Queryname sort the unmapped BAM in parallel with mapped BAM sort.
    SAMTOOLS_QNAME_SORT_UNMAPPED(
        FGUMI_SNAP_ZIPPER_RUN.out.unmapped_bam
            .join(
                ch_meta_unmapped_index_fasta_dict.map { meta, _unmapped_bam, _index, fasta, _dict -> [meta, fasta] },
                by: 0,
            )
            .map { meta, unmapped_bam, fasta -> [meta, unmapped_bam, fasta] },
        null
    )

    // Sort mapped alignments by queryname and emit SAM for zipper stdin.
    SAMTOOLS_QNAME_SORT_MAPPED(
        FGUMI_SNAP_ZIPPER_RUN.out.mapped_bam
            .join(
                ch_meta_unmapped_index_fasta_dict.map { meta, _unmapped_bam, _index, fasta, _dict -> [meta, fasta] },
                by: 0,
            )
            .map { meta, mapped_bam, fasta -> [meta, mapped_bam, fasta] },
        null
    )

    FGUMI_ZIPPER(
        SAMTOOLS_QNAME_SORT_MAPPED.out.sam
            .join(SAMTOOLS_QNAME_SORT_UNMAPPED.out.bam, by: 0)
            .join(
                ch_meta_unmapped_index_fasta_dict.map { meta, _unmapped_bam, _index, fasta, dict -> [meta, fasta, dict] },
                by: 0,
            )
            .map { meta, mapped_sam, unmapped_qname_bam, fasta, dict -> [meta, mapped_sam, unmapped_qname_bam, fasta, dict] }
    )

    FGUMI_TEMPLATE_SORT(FGUMI_ZIPPER.out.bam)

    emit:
    bam            = FGUMI_TEMPLATE_SORT.out.bam
    versions_fgumi = FGUMI_SNAP_ZIPPER_RUN.out.versions_fgumi.mix(FGUMI_ZIPPER.out.versions_fgumi).mix(FGUMI_TEMPLATE_SORT.out.versions_fgumi)
}
