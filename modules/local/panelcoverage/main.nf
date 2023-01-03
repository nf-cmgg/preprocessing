process PANELCOVERAGE {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::pandas=1.4.3 bioconda::bedtools=2.30.0"
    container "cmgg/panelcoverage:latest"

    input:
    tuple val(meta), path(per_base_bed)
    path panel_bed_dir

    output:
    tuple val(meta), path('*.intersected.bed')      , emit: regions_bed
    tuple val(meta), path('*.region.dist.txt')      , emit: regions_txt
    path  "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    for f in ${panel_bed_dir}/*.bed; do
        bedtools intersect \\
            -a ${per_base_bed} \\
            -b \$f \\
            > ${prefix}.intersected.bed

        create_cumulative_region_dist.py \\
            --bed ${prefix}.intersected.bed \\
            > ${prefix}.region.dist.txt
    done


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        pandas: \$(python -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """
}
