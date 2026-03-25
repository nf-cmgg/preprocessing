process SNAPALIGNER_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snap-aligner:2.0.5--h077b44d_2':
        'biocontainers/snap-aligner:2.0.5--h077b44d_2' }"

    input:
    tuple val(meta) , path(reads, stageAs: "?/*"), path(index)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.bai"), optional: true, emit: bai
    tuple val("${task.process}"), val('snap-aligner'), eval("snap-aligner 2>&1 | sed 's/^.*version //;s/.\$//;q'"), topic: versions, emit: versions_snapaligner

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def subcmd = meta.single_end ? "single" : "paired"
    def has_interleaved_flag = args.contains('-pairedInterleavedFastq')
    def args_before_reads = has_interleaved_flag ? '-pairedInterleavedFastq' : ''
    def args_after_reads = has_interleaved_flag ? args.replace('-pairedInterleavedFastq', '').trim() : args

    """
    INDEX=`dirname \$(find -L ./ -name "OverflowTable*")`
    [ -z "\$INDEX" ] && echo "Snap index files not found" 1>&2 && exit 1

    snap-aligner ${subcmd} \
        \$INDEX \
        ${args_before_reads} \
        ${reads} \
        -o ${prefix}.bam \
        -t ${task.cpus} \
        ${args_after_reads}
    """

    stub:
    """
    touch test.bam
    touch test.bam.bai
    """
}
