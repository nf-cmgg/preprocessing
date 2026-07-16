process FGUMI_SNAPZIPSORT {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/05/057cc55ab35ff976996184621e292608f3f53e5fae66aeb36a6ec13659ad6beb/data'
        : 'community.wave.seqera.io/library/fgumi_snap-aligner:fa44bec655a3a203'}"

    input:
    tuple val(meta), path(unmapped_bam), path(index, stageAs: "index/*"), path(fasta), path(fai), path(dict)

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam
    tuple val(meta), path("${prefix}.bam.bai"), emit: bai, optional: true
    tuple val("${task.process}"), val('fgumi'), eval("fgumi --version | sed 's/^fgumi //;q'"), topic: versions
    tuple val("${task.process}"), val('snap-aligner'), eval("snap-aligner 2>&1 | sed 's/^.*version //;s/.\$//;q'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def args4 = task.ext.args4 ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.fgumi"

    """
    # SNAP index directory is resolved from staged index content.
    INDEX_FILE=\$(find -L ./ -name "OverflowTable*" -print -quit)
    [ -z "\$INDEX_FILE" ] && echo "Snap index files not found" 1>&2 && exit 1
    INDEX=\$(dirname "\$INDEX_FILE")

    snap-aligner paired \\
        \$INDEX \\
        ${unmapped_bam} \\
        -t ${task.cpus} \\
        -o -bam - \\
        ${args} \\
    | fgumi sort \\
        --input - \\
        --output ${prefix}.intermediate.bam \\
        --threads ${task.cpus} \\
        ${args2}

    fgumi zipper \\
        --input ${prefix}.intermediate.bam \\
        --unmapped ${unmapped_bam} \\
        --reference ${fasta} \\
        --threads ${task.cpus} \\
        ${args3} \\
    | fgumi sort \\
        --input - \\
        --output ${prefix}.bam \\
        --threads ${task.cpus} \\
        ${args4}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.fgumi"
    """
    touch ${prefix}.bam
    """
}
