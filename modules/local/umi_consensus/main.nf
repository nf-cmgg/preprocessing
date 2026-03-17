include { SAMTOOLS_COLLATE as UMI_SAMTOOLS_COLLATE } from '../../../modules/nf-core/samtools/collate/main.nf'
include { SAMTOOLS_FIXMATE as UMI_SAMTOOLS_FIXMATE } from '../../../modules/nf-core/samtools/fixmate/main.nf'
include { SAMTOOLS_SORT as UMI_SAMTOOLS_SORT_TEMPLATE } from '../../../modules/nf-core/samtools/sort/main.nf'
include { SAMTOOLS_FASTQ as UMI_SAMTOOLS_FASTQ } from '../../../modules/nf-core/samtools/fastq/main.nf'
include { SAMTOOLS_SORT as UMI_SAMTOOLS_SORT_FINAL } from '../../../modules/nf-core/samtools/sort/main.nf'

include { FGBIO_COPYUMIFROMREADNAME as UMI_FGBIO_COPYUMIFROMREADNAME } from '../../../modules/nf-core/fgbio/copyumifromreadname/main.nf'
include { FGBIO_GROUPREADSBYUMI as UMI_FGBIO_GROUPREADSBYUMI } from '../../../modules/nf-core/fgbio/groupreadsbyumi/main.nf'
include { FGBIO_CALLMOLECULARCONSENSUSREADS as UMI_FGBIO_CALLMOLECULARCONSENSUSREADS } from '../../../modules/nf-core/fgbio/callmolecularconsensusreads/main.nf'
include { FGBIO_FILTERCONSENSUSREADS as UMI_FGBIO_FILTERCONSENSUSREADS } from '../../../modules/nf-core/fgbio/filterconsensusreads/main.nf'
include { FGBIO_ZIPPERBAMS as UMI_FGBIO_ZIPPERBAMS } from '../../../modules/nf-core/fgbio/zipperbams/main.nf'

// Local lightweight filtering step used before module-based UMI processing.
// Keeps primary mapped reads and removes unmapped/secondary/supplementary records.
process UMI_LOCAL_SAMTOOLS_VIEW {
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
