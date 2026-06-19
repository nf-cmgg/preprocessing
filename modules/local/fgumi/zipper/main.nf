process FGUMI_ZIPPER {
    tag "$meta.id"
    label 'process_high'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/763a833519c23555be888065f492215f57344155106972e272a0f8df78c57659/data'
        : 'community.wave.seqera.io/library/fgumi_samtools_snap-aligner:c985f9394623a414'}"

    input:
    tuple val(meta), path(mapped_sam), path(unmapped_qname_bam), path(fasta), path(dict)

    output:
    tuple val(meta), path("${prefix}.zipper.bam"), emit: bam
    tuple val("${task.process}"), val('fgumi'), eval("fgumi --version | sed 's/^fgumi //;q'"), topic: versions, emit: versions_fgumi

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.fgumi"

    """
    # mapped_sam and unmapped_qname_bam must be queryname-sorted in the same order.
    fgumi zipper \
            --unmapped ${unmapped_qname_bam} \
            --reference ${fasta} \
            ${args} \
            --output ${prefix}.zipper.bam \
        < ${mapped_sam}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.fgumi"
    """
    touch ${prefix}.zipper.bam
    """
}
