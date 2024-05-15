process PANELCOVERAGE {
    tag "$meta.id - $genelist"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_1' :
        'biocontainers/bedtools:2.31.1--hf5e1c6e_1' }"

    input:
    tuple val(meta), path(perbase), path(perbase_index)
    each path(genelist)

    output:
    tuple val(meta), path("*.mosdepth.region.dist.txt"), emit: regiondist
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cmgg_genelists regiondist --samplename ${prefix} --perbase ${perbase} --genelist ${genelist}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cmgg_genelists: \$(cmgg_genelists -v 2>&1 | sed 's/^.*cmgg_genelists version //')
        bedtools: \$(echo \$(bedtools --version 2>&1) | sed 's/^.*bedtools v//' ))
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_genelist.mosdepth.region.dist.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cmgg_genelists: \$(cmgg_genelists --version 2>&1 | sed 's/^.*cmgg_genelists version //')
        bedtools: \$(echo \$(bedtools --version 2>&1) | sed 's/^.*bedtools v//' ))
    END_VERSIONS
    """
}
