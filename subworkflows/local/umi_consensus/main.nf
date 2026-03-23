include { FASTQ_ALIGN_DNA } from '../../nf-core/fastq_align_dna/main'
include { FASTQ_ALIGN_DNA as FASTQ_ALIGN_DNA_CONSENSUS } from '../../nf-core/fastq_align_dna/main'

include { SAMTOOLS_FASTQ as UMI_SAMTOOLS_FASTQ } from '../../../modules/nf-core/samtools/fastq/main.nf'
include { SAMTOOLS_SORT as UMI_SAMTOOLS_SORT_FINAL } from '../../../modules/nf-core/samtools/sort/main.nf'

include { FGBIO_COPYUMIFROMREADNAME as UMI_FGBIO_COPYUMIFROMREADNAME } from '../../../modules/nf-core/fgbio/copyumifromreadname/main.nf'
include { FGBIO_GROUPREADSBYUMI as UMI_FGBIO_GROUPREADSBYUMI } from '../../../modules/nf-core/fgbio/groupreadsbyumi/main.nf'
include { FGBIO_CALLMOLECULARCONSENSUSREADS as UMI_FGBIO_CALLMOLECULARCONSENSUSREADS } from '../../../modules/nf-core/fgbio/callmolecularconsensusreads/main.nf'
include { FGBIO_FILTERCONSENSUSREADS as UMI_FGBIO_FILTERCONSENSUSREADS } from '../../../modules/nf-core/fgbio/filterconsensusreads/main.nf'
include { FGBIO_ZIPPERBAMS as UMI_FGBIO_ZIPPERBAMS } from '../../../modules/nf-core/fgbio/zipperbams/main.nf'

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

// UMI consensus workflow for DNA samples.
// Input channel shape:
//   [meta, reads, aligner, index, fasta]
// Output channels:
//   bam_bai      -> [meta, bam, bai]
//   family_sizes -> [meta, histogram]
workflow UMI_CONSENSUS_KAPA {
    take:
    ch_meta_reads_aligner_index_fasta // [meta, reads, aligner, index, fasta]

    main:
    // 1) Initial mapping of raw reads with the configured aligner/index.
    FASTQ_ALIGN_DNA(ch_meta_reads_aligner_index_fasta, false)

    // 2) Build reference helper channels reused by downstream modules.
    ch_meta_fasta_fai = ch_meta_reads_aligner_index_fasta
        .map { meta, _reads, _aligner, _index, fasta ->
            def fai = meta.genome_data?.fai ?: '/etc/passwd'
            [meta, fasta, file(fai, checkIfExists: true)]
        }

    ch_meta_fasta = ch_meta_reads_aligner_index_fasta
        .map { meta, _reads, _aligner, _index, fasta -> [meta, fasta] }

    ch_meta_dict = ch_meta_reads_aligner_index_fasta
        .map { meta, _reads, _aligner, _index, _fasta ->
            def dict = meta.genome_data?.dict ?: '/dev/null'
            [meta, file(dict, checkIfExists: true)]
        }

    // 3) Prepare read-pair metadata and UMI tags before consensus calling.
    UMI_SAMTOOLS_PREP_TEMPLATE(
        FASTQ_ALIGN_DNA.out.bam
            .join(FASTQ_ALIGN_DNA.out.bam_index, by: 0)
            .map { meta, bam, bam_index -> [meta, bam, bam_index] },
        ch_meta_fasta_fai
    )

    ch_template_bam_bai = UMI_SAMTOOLS_PREP_TEMPLATE.out.bam_bai

    // 4) Copy UMI from read names to RX tag, then group by UMI families.
    UMI_FGBIO_COPYUMIFROMREADNAME(ch_template_bam_bai)

    UMI_FGBIO_GROUPREADSBYUMI(
        UMI_FGBIO_COPYUMIFROMREADNAME.out.bam,
        channel.value('Adjacency')
    )

    // 5) Call and filter molecular consensus reads.
    UMI_FGBIO_CALLMOLECULARCONSENSUSREADS(
        UMI_FGBIO_GROUPREADSBYUMI.out.bam,
        UMI_FGBIO_GROUPREADSBYUMI.out.bam.map { meta, _bam -> meta.umi_min_reads ?: 2 },
        channel.value(20)
    )

    UMI_FGBIO_FILTERCONSENSUSREADS(
        UMI_FGBIO_CALLMOLECULARCONSENSUSREADS.out.bam,
        ch_meta_fasta,
        UMI_FGBIO_CALLMOLECULARCONSENSUSREADS.out.bam.map { meta, _bam -> meta.umi_min_reads ?: 2 },
        channel.value(45),
        channel.value(0.2)
    )

    UMI_SAMTOOLS_FASTQ(
        UMI_FGBIO_FILTERCONSENSUSREADS.out.bam,
        true
    )

    // 6) Re-map consensus reads with the same aligner/index used upstream.
    FASTQ_ALIGN_DNA_CONSENSUS(
        UMI_SAMTOOLS_FASTQ.out.interleaved
            .join(ch_meta_reads_aligner_index_fasta.map { meta, _reads, aligner, index, fasta -> [meta, aligner, index, fasta] }, by: 0)
            .map { meta, interleaved_fastq, aligner, index, fasta -> [meta, [interleaved_fastq], aligner, index, fasta] }
    , false)

    FASTQ_ALIGN_DNA_CONSENSUS.out.bam
        .join(UMI_FGBIO_FILTERCONSENSUSREADS.out.bam, by: 0)
        .map { meta, mapped_bam, unmapped_bam -> [meta, mapped_bam, unmapped_bam] }
        .set { ch_zipper_bams }

    ch_meta_fasta_fai
        .join(ch_meta_dict, by: 0)
        .map { meta, fasta, fai, dict -> [meta, fasta, fai, dict] }
        .set { ch_zipper_ref }

    // 7) Transfer unmapped metadata back to mapped consensus alignments.
    UMI_FGBIO_ZIPPERBAMS(
        ch_zipper_bams,
        ch_zipper_ref
    )

    // 8) Final coordinate sort + index for downstream CRAM conversion.
    UMI_SAMTOOLS_SORT_FINAL(
        UMI_FGBIO_ZIPPERBAMS.out.bam
            .join(ch_meta_fasta, by: 0)
            .map { meta, bam, fasta -> [meta, bam, fasta] },
        'bai'
    )

    emit:
    bam_bai = UMI_SAMTOOLS_SORT_FINAL.out.bam
        .join(UMI_SAMTOOLS_SORT_FINAL.out.bai, by: 0)
    family_sizes = UMI_FGBIO_GROUPREADSBYUMI.out.histogram
}
