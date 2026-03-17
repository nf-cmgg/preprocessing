include { FASTQ_ALIGN_DNA } from '../../../subworkflows/nf-core/fastq_align_dna/main'
include { FASTQ_ALIGN_DNA as FASTQ_ALIGN_DNA_CONSENSUS } from '../../../subworkflows/nf-core/fastq_align_dna/main'

process UMI_SAMTOOLS_VIEW {
    tag "${meta.id}"
    label 'process_medium'

    container "${task.ext.container ?: (workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0'
        : 'quay.io/biocontainers/samtools:1.22.1--h96c455f_0')}"

    input:
    tuple val(meta), path(sam)

    output:
    tuple val(meta), path("*.umi.filtered.bam"), emit: bam

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools view -F 260 -bh ${sam} > ${prefix}.umi.filtered.bam
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.umi.filtered.bam
    """
}

process UMI_SAMTOOLS_COLLATE {
    tag "${meta.id}"
    label 'process_medium'

    container "${task.ext.container ?: (workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0'
        : 'quay.io/biocontainers/samtools:1.22.1--h96c455f_0')}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.umi.collated.bam"), emit: bam

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools collate -@ ${task.cpus} -O -u -T ${prefix}.collate ${bam} > ${prefix}.umi.collated.bam
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.umi.collated.bam
    """
}

process UMI_SAMTOOLS_FIXMATE {
    tag "${meta.id}"
    label 'process_medium'

    container "${task.ext.container ?: (workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0'
        : 'quay.io/biocontainers/samtools:1.22.1--h96c455f_0')}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.umi.fixmate.bam"), emit: bam

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools fixmate -m -@ ${task.cpus} ${bam} ${prefix}.umi.fixmate.bam
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.umi.fixmate.bam
    """
}

process UMI_SAMTOOLS_TEMPLATE_SORT {
    tag "${meta.id}"
    label 'process_medium'

    container "${task.ext.container ?: (workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0'
        : 'quay.io/biocontainers/samtools:1.22.1--h96c455f_0')}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.umi.template_sorted.bam"), emit: bam

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools sort --template-coordinate -@ ${task.cpus} -T ${prefix}.sort -o ${prefix}.umi.template_sorted.bam ${bam}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.umi.template_sorted.bam
    """
}

process UMI_FGBIO_COPY_UMI {
    tag "${meta.id}"
    label 'process_medium'

    container "${task.ext.container ?: (workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/fgbio:2.1.0--hdfd78af_0'
        : 'quay.io/biocontainers/fgbio:2.1.0--hdfd78af_0')}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.umi_mapped_filtered.bam"), emit: bam

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3072
    if (!task.memory) {
        log.info('[UMI_FGBIO_COPY_UMI] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    JAVA_TOOL_OPTIONS="-Xmx${avail_mem}M" fgbio --tmp-dir "${task.workDir}" --compression 1 --async-io CopyUmiFromReadName --input ${bam} --output ${prefix}.umi_mapped_filtered.bam
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.umi_mapped_filtered.bam
    """
}

process UMI_FGBIO_GROUP_READS {
    tag "${meta.id}"
    label 'process_medium'

    container "${task.ext.container ?: (workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/fgbio:2.1.0--hdfd78af_0'
        : 'quay.io/biocontainers/fgbio:2.1.0--hdfd78af_0')}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.umi_mapped.grouped.bam"), emit: grouped_bam
    tuple val(meta), path("*.umi_tag-family-sizes_counts.txt"), emit: family_sizes

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3072
    if (!task.memory) {
        log.info('[UMI_FGBIO_GROUP_READS] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    JAVA_TOOL_OPTIONS="-Xmx${avail_mem}M" fgbio --tmp-dir "${task.workDir}" --compression 1 --async-io GroupReadsByUmi --input ${bam} --strategy adjacency --edits 1 -t RX --output ${prefix}.umi_mapped.grouped.bam --family-size-histogram ${prefix}.umi_tag-family-sizes_counts.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.umi_mapped.grouped.bam
    touch ${prefix}.umi_tag-family-sizes_counts.txt
    """
}

process UMI_FGBIO_CALL_CONSENSUS {
    tag "${meta.id}"
    label 'process_medium'

    container "${task.ext.container ?: (workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/fgbio:2.1.0--hdfd78af_0'
        : 'quay.io/biocontainers/fgbio:2.1.0--hdfd78af_0')}"

    input:
    tuple val(meta), path(grouped_bam)

    output:
    tuple val(meta), path("*.umi_consensus.minreads*.raw.unmapped.bam"), emit: raw_unmapped

    script:
    def minReads = meta.umi_min_reads ?: 2
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3072
    if (!task.memory) {
        log.info('[UMI_FGBIO_CALL_CONSENSUS] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    JAVA_TOOL_OPTIONS="-Xmx${avail_mem}M" fgbio --tmp-dir "${task.workDir}" --compression 0 CallMolecularConsensusReads --input ${grouped_bam} --output ${prefix}.umi_consensus.minreads${minReads}.raw.unmapped.bam --error-rate-pre-umi 45 --error-rate-post-umi 40 --min-input-base-quality 20 --min-reads ${minReads} --max-reads 50 --output-per-base-tags false --read-name-prefix consensus --threads ${task.cpus}
    """

    stub:
    def minReads = meta.umi_min_reads ?: 2
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.umi_consensus.minreads${minReads}.raw.unmapped.bam
    """
}

process UMI_FGBIO_FILTER_CONSENSUS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fe/fe9479adc5e6e0a1c125d346fdfa0dd313834249e9c55c40e8d44ec3a48c6559/data' :
        'community.wave.seqera.io/library/fgbio:3.1.1--6c9a88faf1d62b6c' }"

    input:
    tuple val(meta), path(raw_unmapped), path(reference)

    output:
    tuple val(meta), path("*.umi_consensus.minreads*.unmapped.bam"), emit: unmapped

    script:
    def minReads = meta.umi_min_reads ?: 2
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3072
    if (!task.memory) {
        log.info('[UMI_FGBIO_FILTER_CONSENSUS] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    JAVA_TOOL_OPTIONS="-Xmx200g" fgbio --tmp-dir "${task.workDir}" --compression 1 FilterConsensusReads --input ${raw_unmapped} --output ${prefix}.umi_consensus.minreads${minReads}.unmapped.bam --ref ${reference} --min-reads ${minReads} --max-base-error-rate 0.2 --min-base-quality 45
    """

    stub:
    def minReads = meta.umi_min_reads ?: 2
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.umi_consensus.minreads${minReads}.unmapped.bam
    """
}

process UMI_SAMTOOLS_FASTQ {
    tag "${meta.id}"
    label 'process_medium'

    container "${task.ext.container ?: (workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0'
        : 'quay.io/biocontainers/samtools:1.22.1--h96c455f_0')}"

    input:
    tuple val(meta), path(unmapped)

    output:
    tuple val(meta), path("*.umi_consensus.interleaved.fastq"), path(unmapped), emit: fastq_unmapped

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools fastq ${unmapped} > ${prefix}.umi_consensus.interleaved.fastq
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.umi_consensus.interleaved.fastq
    """
}

process UMI_FGBIO_ZIPPER_BAMS {
    tag "${meta.id}"
    label 'process_medium'

    container "${task.ext.container ?: (workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/fgbio:2.1.0--hdfd78af_0'
        : 'quay.io/biocontainers/fgbio:2.1.0--hdfd78af_0')}"

    input:
    tuple val(meta), path(remap_bam), path(unmapped), path(reference)

    output:
    tuple val(meta), path("*.umi_consensus.zippered.bam"), emit: bam

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3072
    if (!task.memory) {
        log.info('[UMI_FGBIO_ZIPPER_BAMS] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    JAVA_TOOL_OPTIONS="-Xmx${avail_mem}M" fgbio --tmp-dir "${task.workDir}" --compression 0 --async-io ZipperBams --input ${remap_bam} --unmapped ${unmapped} --ref ${reference} --tags-to-reverse Consensus --tags-to-revcomp Consensus --output ${prefix}.umi_consensus.zippered.bam
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.umi_consensus.zippered.bam
    """
}

process UMI_SAMTOOLS_FINAL_SORT {
    tag "${meta.id}"
    label 'process_medium'

    container "${task.ext.container ?: (workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0'
        : 'quay.io/biocontainers/samtools:1.22.1--h96c455f_0')}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.umi_consensus.minreads*.mapped.filtered.bam"), path("*.umi_consensus.minreads*.mapped.filtered.bam.bai"), emit: bam_bai
    tuple val("${task.process}"), val('samtools'), eval('samtools --version | sed -n "1p" | sed "s/^samtools //"'), emit: versions_samtools, topic: versions

    script:
    def minReads = meta.umi_min_reads ?: 2
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools sort --threads ${task.cpus} -l 9 -T ${prefix}.consensus_sort -o ${prefix}.umi_consensus.minreads${minReads}.mapped.filtered.bam --write-index ${bam}
    """

    stub:
    def minReads = meta.umi_min_reads ?: 2
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.umi_consensus.minreads${minReads}.mapped.filtered.bam
    touch ${prefix}.umi_consensus.minreads${minReads}.mapped.filtered.bam.bai
    """
}

workflow UMI_CONSENSUS_KAPA {
    take:
    ch_meta_reads_aligner_index_fasta // [meta, reads, aligner, index, fasta]

    main:
    FASTQ_ALIGN_DNA(ch_meta_reads_aligner_index_fasta, false)
    UMI_SAMTOOLS_VIEW(FASTQ_ALIGN_DNA.out.bam)
    UMI_SAMTOOLS_COLLATE(UMI_SAMTOOLS_VIEW.out.bam)
    UMI_SAMTOOLS_FIXMATE(UMI_SAMTOOLS_COLLATE.out.bam)
    UMI_SAMTOOLS_TEMPLATE_SORT(UMI_SAMTOOLS_FIXMATE.out.bam)
    UMI_FGBIO_COPY_UMI(UMI_SAMTOOLS_TEMPLATE_SORT.out.bam)
    UMI_FGBIO_GROUP_READS(UMI_FGBIO_COPY_UMI.out.bam)

    UMI_FGBIO_CALL_CONSENSUS(UMI_FGBIO_GROUP_READS.out.grouped_bam)

    UMI_FGBIO_FILTER_CONSENSUS(
        UMI_FGBIO_CALL_CONSENSUS.out.raw_unmapped.join(
            ch_meta_reads_aligner_index_fasta.map { meta, _reads, _aligner, _index, fasta -> [meta, fasta] },
            by: 0
        ).map { meta, raw_unmapped, reference -> [meta, raw_unmapped, reference] }
    )

    UMI_SAMTOOLS_FASTQ(UMI_FGBIO_FILTER_CONSENSUS.out.unmapped)

    FASTQ_ALIGN_DNA_CONSENSUS(
        UMI_SAMTOOLS_FASTQ.out.fastq_unmapped
            .join(ch_meta_reads_aligner_index_fasta.map { meta, _reads, aligner, index, fasta -> [meta, aligner, index, fasta] }, by: 0)
            .map { meta, interleaved_fastq, _unmapped, aligner, index, fasta -> [meta, [interleaved_fastq], aligner, index, fasta] }
    , false)

    UMI_FGBIO_ZIPPER_BAMS(
        FASTQ_ALIGN_DNA_CONSENSUS.out.bam
            .join(UMI_SAMTOOLS_FASTQ.out.fastq_unmapped.map { meta, _interleaved_fastq, unmapped -> [meta, unmapped] }, by: 0)
            .join(ch_meta_reads_aligner_index_fasta.map { meta, _reads, _aligner, _index, fasta -> [meta, fasta] }, by: 0)
            .map { meta, remap_bam, unmapped, reference -> [meta, remap_bam, unmapped, reference] }
    )

    UMI_SAMTOOLS_FINAL_SORT(UMI_FGBIO_ZIPPER_BAMS.out.bam)

    emit:
    bam_bai = UMI_SAMTOOLS_FINAL_SORT.out.bam_bai
    family_sizes = UMI_FGBIO_GROUP_READS.out.family_sizes
}
