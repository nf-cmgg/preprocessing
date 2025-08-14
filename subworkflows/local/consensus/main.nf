include { FGBIO_COPYUMIFROMREADNAME } from '../../../modules/nf-core/fgbio/copyumifromreadname/main'
include { FGBIO_CALLDUPLEXCONSENSUSREADS } from '../../../modules/nf-core/fgbio/callduplexconsensusreads/main'
include { FGBIO_CALLMOLECULARCONSENSUSREADS } from '../../../modules/nf-core/fgbio/callmolecularconsensusreads/main'
include { FGBIO_COLLECTDUPLEXSEQMETRICS } from '../../../modules/nf-core/fgbio/collectduplexseqmetrics/main'
include { FGBIO_FASTQTOBAM } from '../../../modules/nf-core/fgbio/fastqtobam/main'
include { FGBIO_FILTERCONSENSUSREADS } from '../../../modules/nf-core/fgbio/filterconsensusreads/main'
include { FGBIO_GROUPREADSBYUMI } from '../../../modules/nf-core/fgbio/groupreadsbyumi/main'
include { FGBIO_SORTBAM } from '../../../modules/nf-core/fgbio/sortbam/main'
include { FGBIO_ZIPPERBAMS } from '../../../modules/nf-core/fgbio/zipperbams/main'

workflow CONSENSUS {
    take:
    ch_fastq                // channel: tuple(meta, fastq1, fastq2) for SE/PE/duplex samples
    ch_reference            // channel: reference genome fasta file
    ch_umi_in_readname      // boolean

    main:
    ch_fastq
    .combine(ch_umi_in_readname)
    .branch { tuple ->
        meta, fastq1, fastq2, umi_in_readname = tuple
        umi_in_readname: umi_in_readname == true
        umi_in_sequence: umi_in_readname == false
    }
    .set { ch_fastq_branch }

    // UMI is in read name: run CopyUmiFromReadName before FastqToBam
    ch_fastq_branch.umi_in_readname
    .map { meta, fastq1, fastq2 ->
        // Step 1: Extract UMI from read name and copy to BAM tag
        COPYUMIFROMREADNAME([meta, fastq1, fastq2])
        // Step 2: Convert FASTQ to BAM (UMI tag already present)
        FASTQTOBAM(COPYUMIFROMREADNAME.out.bam, ch_reference, ch_umi_info)
        // Step 3: Sort BAM
        SORTBAM(FASTQTOBAM.out.bam)
        // Step 4: Group reads by UMI
        GROUPREADSBYUMI(SORTBAM.out.bam)
        // Step 5: Call molecular consensus reads
        CALLMOLECULARCONSENSUSREADS(GROUPREADSBYUMI.out.bam)
        // Step 6: If duplex UMI, call duplex consensus reads
        if (meta.umi_kit == 'DUPLEX') {
            CALLDUPLEXCONSENSUSREADS(CALLMOLECULARCONSENSUSREADS.out.bam)
        }
        // Step 7: Filter consensus reads
        FILTERCONSENSUSREADS(...)
        // Additional steps as needed
    }

    // UMI is in sequence: run FastqToBam directly
    ch_fastq_branch.umi_in_sequence
    .map { meta, fastq1, fastq2 ->
        // Step 1: Convert FASTQ to BAM (UMI extracted from sequence/index)
        FASTQTOBAM([meta, fastq1, fastq2], ch_reference, ch_umi_info)
        // Step 2: Sort BAM
        SORTBAM(FASTQTOBAM.out.bam)
        // Step 3: Group reads by UMI
        GROUPREADSBYUMI(SORTBAM.out.bam)
        // Step 4: Call molecular consensus reads
        CALLMOLECULARCONSENSUSREADS(GROUPREADSBYUMI.out.bam)
        // Step 5: If duplex UMI, call duplex consensus reads
        if (meta.umi_kit == 'DUPLEX') {
            CALLDUPLEXCONSENSUSREADS(CALLMOLECULARCONSENSUSREADS.out.bam)
        }
        // Step 6: Filter consensus reads
        FILTERCONSENSUSREADS(...)
        // Additional steps as needed
    }

    emit:
    consensus_bam = ... // Final consensus BAM/CRAM output
}
