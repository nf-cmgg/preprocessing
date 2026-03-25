process UMI_SAMTOOLS_STRIP_PG {
    tag "$meta.id"
    label 'process_single'

    conda "${projectDir}/modules/nf-core/samtools/view/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0' :
        'biocontainers/samtools:1.22.1--h96c455f_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.no_pg"

    """
    samtools view -H ${bam} | awk '\$1 != "@PG"' > ${prefix}.header.sam
    samtools reheader --no-PG ${prefix}.header.sam ${bam} > ${prefix}.bam
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.no_pg"
    """
    touch ${prefix}.bam
    """
}
