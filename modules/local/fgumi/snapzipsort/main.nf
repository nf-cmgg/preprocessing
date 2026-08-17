process FGUMI_SNAPZIPSORT {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d9/d9ad0d2e6f34163511be159dcdb92b989fb6cd3f3eebe5bc0c943a3318503a96/data'
        : 'community.wave.seqera.io/library/fgumi_snap-aligner:1dd474f4076b25af'}"

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

    fgumi fastq \\
        --input ${unmapped_bam} \\
        ${args} \\
    | snap-aligner paired \\
        \$INDEX \\
        -pairedInterleavedFastq - \\
        -t ${task.cpus} \\
        -o -bam - \\
        ${args2} \\
    | fgumi zipper \\
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
