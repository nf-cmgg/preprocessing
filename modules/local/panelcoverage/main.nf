process PANELCOVERAGE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_1'
        : 'biocontainers/bedtools:2.31.1--hf5e1c6e_1'}"

    input:
    tuple val(meta), path(perbase), path(perbase_index), path(genelists)

    output:
    tuple val(meta), path("*.mosdepth.region.dist.txt"), emit: regiondist
    tuple val("${task.process}"), val('cmgg_genelists'), eval('cmgg_genelists -v 2>&1 | sed \"s/^.*cmgg_genelists version //\"'), emit: versions_cmgg_genelists, topic: versions
    tuple val("${task.process}"), val('bedtools'), eval('bedtools --version 2>&1 | sed \"s/^.*bedtools v//\"'), emit: versions_bedtools, topic: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${genelists} | tr ' ' '\n' | xargs -n 1 -P ${task.cpus} -I {} cmgg_genelists regiondist --samplename ${prefix} --perbase ${perbase} --genelist {}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    for GENELIST in ${genelists}
    do
        name=\$(basename \$GENELIST .bed)
        touch ${prefix}_\${name}.mosdepth.region.dist.txt
    done
    """
}
