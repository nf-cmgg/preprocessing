process FGUMI_FILTER {
    tag "$meta.id"
    label 'process_medium'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/95/954170443a820787c9e02ef2135ebb8ec29c6b03633b0d61b5fafa98c59a1cce/data'
        : 'community.wave.seqera.io/library/fgumi_r-base_r-ggplot2_r-scales:09c99070b82c1c28'}"

    input:
    tuple val(meta), path(bam), path(fasta)

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam
    tuple val(meta), path("${prefix}.bam.bai"), emit: bai
    tuple val(meta), path("${prefix}.filtering_metrics.txt"), optional: true, emit: filtering_metrics
    tuple val("${task.process}"), val('fgumi'), eval("fgumi --version | sed 's/^fgumi //;q'"), topic: versions, emit: versions_fgumi

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def sort_args = task.ext.args2 ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.fgumi.filter"

    """
    fgumi filter \
        --input ${bam} \
        --output ${prefix}.filtered.bam \
        --ref ${fasta} \
        ${args}

    fgumi sort \
        --input ${prefix}.filtered.bam \
        --output ${prefix}.bam \
        --order coordinate \
        --write-index \
        ${sort_args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.fgumi.filter"
    """
    touch ${prefix}.bam
    touch ${prefix}.bam.bai
    touch ${prefix}.filtering_metrics.txt
    """
}
