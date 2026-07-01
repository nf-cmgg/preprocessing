process FGUMI_EXTRACT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0e/0efc3d689f6d1cbfb55ae9a20ce8758e631c6cc6aca51186f9e82568542a3896/data' :
          'community.wave.seqera.io/library/fgumi:0.4.0--9fce7cd93e5aceee' }"

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
