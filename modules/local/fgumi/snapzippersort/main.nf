process FGUMI_SNAP_ZIPPER_SORT {
    tag "$meta.id"
    label 'process_high'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/24/2466270633749543330f352e34588f142de4988585ce63e7f22ee5ed1ff57450/data'
        : 'community.wave.seqera.io/library/fgumi_samtools_snap-aligner:c9ba911435350668'}"

    input:
    tuple val(meta), path(unmapped_bam), path(index, stageAs: "index/*"), path(fasta), path(dict)

    output:
    tuple val(meta), path("${prefix}.template.bam"), emit: bam
    tuple val("${task.process}"), val('fgumi'), eval("fgumi --version | sed 's/^fgumi //;q'"), topic: versions, emit: versions_fgumi

    when:
    task.ext.when == null || task.ext.when

    script:
    def snap_args = task.ext.args ?: ''
    def zipper_args = task.ext.args2 ?: ''
    def sort_args = task.ext.args3 ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.fgumi"

    """
    INDEX=`dirname \$(find -L ./ -name "OverflowTable*" | head -n1)`
    [ -z "\$INDEX" ] && echo "Snap index files not found" 1>&2 && exit 1

    # Ensure zipper and fastq read exactly the same queryname-ordered unmapped stream.
    samtools sort \
        -n \
        -@ ${task.cpus} \
        -m 1G \
        -o ${prefix}.unmapped.queryname.bam \
        ${unmapped_bam}

    fgumi fastq --input ${prefix}.unmapped.queryname.bam \
        | snap-aligner paired \
            \$INDEX \
            -pairedInterleavedFastq - \
            -o -sam - \
            -t ${task.cpus} \
            ${snap_args} \
        | samtools sort \
            -n \
            -@ ${task.cpus} \
            -m 1G \
            -O SAM \
            - \
        | fgumi zipper \
            --unmapped ${prefix}.unmapped.queryname.bam \
            --reference ${fasta} \
            ${zipper_args} \
        | fgumi sort \
            --input /dev/stdin \
            --output ${prefix}.template.bam \
            --order template-coordinate \
            ${sort_args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.fgumi"
    """
    touch ${prefix}.template.bam
    """
}
