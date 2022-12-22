process ELPREP_SFM {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::elprep=5.1.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/elprep:5.1.2--he881be0_0':
        'quay.io/biocontainers/elprep:5.1.2--he881be0_0' }"

    input:
    tuple val(meta), path(bam, stageAs: "input/*")

    output:
    tuple val(meta), path("output/**.{bam,sam}")    ,emit: bam
    tuple val(meta), path("*.metrics.txt")          ,optional: true, emit: metrics
    path "versions.yml"                             ,emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = args.contains("--output-type sam") ? "sam" : "bam"

    // single end
    def single_end_cmd = meta.single_end ? "--single-end" : ""

    """
    elprep sfm ./input output/${prefix}.${suffix} \\
        ${single_end_cmd} \\
        --nr-of-threads ${task.cpus} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        elprep: \$(elprep 2>&1 | head -n2 | tail -n1 |sed 's/^.*version //;s/ compiled.*\$//')
    END_VERSIONS
    """
}
