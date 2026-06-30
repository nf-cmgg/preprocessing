process FGUMI_ZIPPER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c5/c566f9e20f9eb4c5be9ff5a68e854f974caae916d67b4e03eb30eece186b73e8/data'
        : 'community.wave.seqera.io/library/fgumi_samtools_snap-aligner:1708ad8bd6e764b6'}"

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
