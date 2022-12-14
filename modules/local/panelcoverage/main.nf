process PANELCOVERAGE {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? 'bioconda::pandas=1.4.3' : null)
    container "quay.io/biocontainers/pandas:1.4.3"

    input:
    tuple val(meta), path(per_base_panel_bed)

    output:
    tuple val(meta), path('*.region.dist.txt')      , optional:true, emit: regions_txt
    path  "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    create_cumulative_region_dist.py \\
        --bed ${per_base_panel_bed} \\
        > ${prefix}.region.dist.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pandas: \$(python -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """
}
