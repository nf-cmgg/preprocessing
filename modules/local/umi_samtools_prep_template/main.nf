process UMI_SAMTOOLS_PREP_TEMPLATE {
    tag "$meta.id"
    label 'process_medium'

    conda "${projectDir}/modules/nf-core/samtools/view/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0' :
        'biocontainers/samtools:1.22.1--h96c455f_0' }"

    input:
    tuple val(meta), path(bam), path(bam_index)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("${prefix}.bam"), path("${prefix}.bam.bai"), emit: bam_bai
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def view_args = task.ext.args ?: '-F 260 -bh --output-fmt bam'
    def collate_args = task.ext.args2 ?: '-O -u --output-fmt bam'
    def fixmate_args = task.ext.args3 ?: '-m --output-fmt bam'
    def sort_args = task.ext.args4 ?: '--template-coordinate --output-fmt bam'
    prefix = task.ext.prefix ?: "${meta.id}.template"
    def sort_memory = (task.memory.mega / task.cpus * 0.75).intValue()
    """
    test -f ${bam_index}

    samtools view \
        ${view_args} \
        --threads ${task.cpus} \
        ${bam} \
    | samtools collate \
        ${collate_args} \
        --threads ${task.cpus} \
        - \
    | samtools fixmate \
        ${fixmate_args} \
        --threads ${task.cpus} \
        - - \
    | samtools sort \
        ${sort_args} \
        --threads ${task.cpus} \
        --reference ${fasta} \
        -m ${sort_memory}M \
        -o ${prefix}.bam##idx##${prefix}.bam.bai --write-index \
        -
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.template"
    """
    touch ${prefix}.bam
    touch ${prefix}.bam.bai
    """
}
