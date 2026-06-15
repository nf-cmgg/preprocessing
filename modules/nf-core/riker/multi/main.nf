process RIKER_MULTI {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9b/9b1fa220679fac2077935260a88165214d8625b813ecb39cf129aecfba982c65/data' :
        'community.wave.seqera.io/library/riker:0.2.0--aa377939e8395424' }"

    input:
    tuple val(meta), path(bam), path(bai), path(roi), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.{txt,pdf}"), emit: metrics
    tuple val("${task.process}"), val('riker'), eval("riker --version 2>&1 | sed 's/riker //'") , topic: versions, emit: versions_riker

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def ref         = fasta ? "-r ${fasta}" : ''
    def hybcap_opts = roi ? "--hybcap::baits ${roi} --hybcap::targets ${roi}" : ''
    """
    riker multi \\
        -i ${bam} \\
        ${ref} \\
        -o ${prefix} \\
        --threads ${task.cpus} \\
        ${hybcap_opts} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.alignment-metrics.txt
    touch ${prefix}.base-distribution-by-cycle.txt
    touch ${prefix}.mean-quality-by-cycle.txt
    touch ${prefix}.quality-score-distribution.txt
    touch ${prefix}.error-mismatch.txt
    touch ${prefix}.error-overlap.txt
    touch ${prefix}.error-indel.txt
    touch ${prefix}.gcbias-detail.txt
    touch ${prefix}.gcbias-summary.txt
    touch ${prefix}.hybcap-metrics.txt
    touch ${prefix}.hybcap-per-target.txt
    touch ${prefix}.hybcap-per-base.txt
    touch ${prefix}.isize-metrics.txt
    touch ${prefix}.isize-histogram.txt
    touch ${prefix}.wgs-metrics.txt
    touch ${prefix}.wgs-coverage.txt
    touch ${prefix}.base-distribution-by-cycle.pdf
    touch ${prefix}.gcbias-chart.pdf
    touch ${prefix}.isize-histogram.pdf
    touch ${prefix}.mean-quality-by-cycle.pdf
    touch ${prefix}.quality-score-distribution.pdf
    touch ${prefix}.wgs-coverage.pdf
    """
}
