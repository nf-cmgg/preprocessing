process FGUMI_EXTRACT {
    tag "$meta.id"
    label 'process_medium'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/95/954170443a820787c9e02ef2135ebb8ec29c6b03633b0d61b5fafa98c59a1cce/data'
        : 'community.wave.seqera.io/library/fgumi_r-base_r-ggplot2_r-scales:09c99070b82c1c28'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam
    tuple val("${task.process}"), val('fgumi'), eval("fgumi --version | sed 's/^fgumi //;q'"), topic: versions, emit: versions_fgumi

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.fgumi.unmapped"
    // Keep module portable: only meta.id is assumed, with optional task.ext overrides.
    def sample_name = task.ext.sample_name ?: meta.id
    def library_name = task.ext.library_name ?: meta.id

    """
    fgumi extract \
        --inputs ${reads} \
        --output ${prefix}.bam \
        --sample "${sample_name}" \
        --library "${library_name}" \
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.fgumi.unmapped"
    """
    touch ${prefix}.bam
    """
}
