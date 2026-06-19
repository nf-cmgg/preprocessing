process FGUMI_SNAP_ALIGN {
    tag "$meta.id"
    label 'process_high'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/763a833519c23555be888065f492215f57344155106972e272a0f8df78c57659/data'
        : 'community.wave.seqera.io/library/fgumi_samtools_snap-aligner:c985f9394623a414'}"

    input:
    tuple val(meta), path(unmapped_bam), path(index, stageAs: "index/*"), path(fasta), path(dict)

    output:
    tuple val(meta), path("${prefix}.snap.bam"), emit: mapped_bam
    tuple val(meta), path(unmapped_bam), emit: unmapped_bam
    tuple val("${task.process}"), val('fgumi'), eval("fgumi --version | sed 's/^fgumi //;q'"), topic: versions, emit: versions_fgumi

    when:
    task.ext.when == null || task.ext.when

    script:
    def snap_args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.fgumi"

    """
    # SNAP index directory is resolved from staged index content.
    INDEX_FILE=\$(find -L ./ -name "OverflowTable*" -print -quit)
    [ -z "\$INDEX_FILE" ] && echo "Snap index files not found" 1>&2 && exit 1
    INDEX=\$(dirname "\$INDEX_FILE")

    fgumi fastq --input ${unmapped_bam} \
        | snap-aligner paired \
            \$INDEX \
            -pairedInterleavedFastq - \
            -o ${prefix}.snap.bam \
            -t ${task.cpus} \
            ${snap_args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.fgumi"
    """
    touch ${prefix}.snap.bam
    touch ${unmapped_bam}
    """
}
