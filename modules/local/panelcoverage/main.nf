process PANELCOVERAGE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_1' :
        'biocontainers/bedtools:2.31.1--hf5e1c6e_1' }"

    input:
    tuple val(meta), path(perbase), path(perbase_index)
    path(genelists)

    output:
    tuple val(meta), path("*.mosdepth.region.dist.txt"), emit: regiondist
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    for GENELIST in $genelists
    do
        cmgg_genelists regiondist --samplename ${prefix} --perbase ${perbase} --genelist \$GENELIST
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cmgg_genelists: \$(cmgg_genelists -v 2>&1 | sed 's/^.*cmgg_genelists version //')
        bedtools: \$(echo \$(bedtools --version 2>&1) | sed 's/^.*bedtools v//' ))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    for GENELIST in $genelists
    do
        name=\$(basename \$GENELIST .bed)
        touch ${prefix}_\${name}.mosdepth.region.dist.txt
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cmgg_genelists: \$(cmgg_genelists --version 2>&1 | sed 's/^.*cmgg_genelists version //')
        bedtools: \$(echo \$(bedtools --version 2>&1) | sed 's/^.*bedtools v//' ))
    END_VERSIONS
    """
}
