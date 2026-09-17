process FGUMI_EXTRACTSORT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/99/99c92db2efcbcc4d20f2541060e09ed0f5a4338927f4a2ae3ab683de585efeb6/data'
        : 'community.wave.seqera.io/library/fgumi:0.7.0--d91f99b4cd4aae5a'}"

    input:
    tuple val(meta), path(reads), val(library)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('fgumi'), eval('fgumi --version | sed "s/^fgumi //"'), topic: versions, emit: versions_fgumi

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    fgumi extract \\
        --inputs ${reads.join(' ')} \\
        --output - \\
        --sample ${prefix} \\
        --library "${library}" \\
        ${args} \\
    | fgumi sort \\
        --input - \\
        --output ${prefix}.bam \\
        --threads ${task.cpus} \\
        --max-memory ${task.memory.toGiga() / task.cpus}G \\
        --tmp-dir . \\
        ${args2}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    """
}
