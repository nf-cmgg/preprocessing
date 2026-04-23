process FGUMI_GROUP {
    tag "$meta.id"
    label 'process_medium'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f7267104b209869695781a5f4585c490a61250269f8c6f14068535a3962b865a/data'
        : 'community.wave.seqera.io/library/fgumi:0.2.0--fe028e7a64e5da27'}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam
    tuple val(meta), path("${prefix}.grouping_metrics.txt"), optional: true, emit: grouping_metrics
    tuple val(meta), path("${prefix}.family_size_histogram.txt"), optional: true, emit: family_size_histogram
    tuple val("${task.process}"), val('fgumi'), eval("fgumi --version | sed 's/^fgumi //;q'"), topic: versions, emit: versions_fgumi

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Derive per-thread queue memory from requested process resources.
    def queue_memory_mb = (task.memory.mega / task.cpus * 0.75).intValue()
    prefix = task.ext.prefix ?: "${meta.id}.fgumi.group"

    """
    fgumi group \
        --input ${bam} \
        --output ${prefix}.bam \
        --threads ${task.cpus} \
        --queue-memory ${queue_memory_mb} \
        --queue-memory-per-thread \
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.fgumi.group"
    """
    touch ${prefix}.bam
    touch ${prefix}.grouping_metrics.txt
    touch ${prefix}.family_size_histogram.txt
    """
}
