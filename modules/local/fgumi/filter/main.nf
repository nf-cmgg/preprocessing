process FGUMI_FILTER {
    tag "$meta.id"
    label 'process_medium'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f7267104b209869695781a5f4585c490a61250269f8c6f14068535a3962b865a/data'
        : 'community.wave.seqera.io/library/fgumi:0.2.0--fe028e7a64e5da27'}"

    input:
    tuple val(meta), path(bam), path(fasta)

    output:
    tuple val(meta), path("${prefix}.filtered.bam"), emit: bam
    tuple val(meta), path("${prefix}.filtering_metrics.txt"), optional: true, emit: filtering_metrics
    tuple val("${task.process}"), val('fgumi'), eval("fgumi --version | sed 's/^fgumi //;q'"), topic: versions, emit: versions_fgumi

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.fgumi.filter"

    """
    fgumi filter \
        --input ${bam} \
        --output ${prefix}.filtered.bam \
        --ref ${fasta} \
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.fgumi.filter"
    """
    touch ${prefix}.filtered.bam
    touch ${prefix}.filtering_metrics.txt
    """
}
