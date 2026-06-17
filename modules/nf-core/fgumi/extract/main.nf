process FGUMI_EXTRACT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/394ef34259ba03c393d25ce2559530fe6df2c6f125a27be69da89ca7a1c70e30/data' :
          'community.wave.seqera.io/library/fgumi:0.3.0--bcbb552cbefcbda1' }"

    input:
    tuple val(meta), path(reads), val(library)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('fgumi'), eval('fgumi --version | sed "s/^fgumi //"'), topic: versions, emit: versions_fgumi

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    fgumi extract \\
        --inputs ${reads.join(' ')} \\
        --output ${prefix}.bam \\
        ${args} \\
        --sample ${prefix} \\
        --library "${library}"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    """
}
