include { FGUMI_SORT as FGUMI_TEMPLATE_SORT } from "../../../nf-core/fgumi/sort/main.nf"

process FGUMI_SNAP_ZIPPER_RUN {
    tag "$meta.id"
    label 'process_high'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/763a833519c23555be888065f492215f57344155106972e272a0f8df78c57659/data'
        : 'community.wave.seqera.io/library/fgumi_samtools_snap-aligner:c985f9394623a414'}"

    input:
    tuple val(meta), path(unmapped_bam), path(index, stageAs: "index/*"), path(fasta), path(dict)

    output:
    tuple val(meta), path("${prefix}.zipper.bam"), emit: bam
    tuple val("${task.process}"), val('fgumi'), eval("fgumi --version | sed 's/^fgumi //;q'"), topic: versions, emit: versions_fgumi

    when:
    task.ext.when == null || task.ext.when

    script:
    def snap_args = task.ext.args ?: ''
    def zipper_args = task.ext.args2 ?: ''
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

    fgumi zipper \
            --unmapped ${unmapped_bam} \
            --reference ${fasta} \
            ${zipper_args} \
            --output ${prefix}.zipper.bam \
        < ${prefix}.snap.sam
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.fgumi"
    """
    touch ${prefix}.zipper.bam
    """
}

workflow FGUMI_SNAP_ZIPPER {
    take:
    ch_meta_unmapped_index_fasta_dict

    main:
    FGUMI_SNAP_ZIPPER_RUN(ch_meta_unmapped_index_fasta_dict)
    FGUMI_TEMPLATE_SORT(FGUMI_SNAP_ZIPPER_RUN.out.bam)

    emit:
    bam            = FGUMI_TEMPLATE_SORT.out.bam
    versions_fgumi = FGUMI_SNAP_ZIPPER_RUN.out.versions_fgumi.mix(FGUMI_TEMPLATE_SORT.out.versions_fgumi)
}
