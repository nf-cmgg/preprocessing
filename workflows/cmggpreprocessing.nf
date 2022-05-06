/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowCmggpreprocessing.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { BAM_STATS_SAMTOOLS} from '../subworkflows/nf-core/subworkflows/bam_stats_samtools'
include { DEMULTIPLEX       } from '../subworkflows/local/demultiplex'
include { INPUT_CHECK       } from '../subworkflows/local/input_check'
inlcude { BAM_QC_PICARD     } from '../subworkflows/nf-core/subworkflows/bam_qc_picard'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'
include { FASTQC                      } from '../modules/nf-core/modules/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/modules/multiqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow CMGGPREPROCESSING {

    ch_versions         = Channel.empty()
    ch_multiqc_files    = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    ch_flowcells        = INPUT_CHECK (ch_input).out.flowcells
    ch_versions         = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // DEMULTIPLEXING
    //
    // SUBWORKFLOW: demultiplex
    ch_raw_fastq        = DEMULTIPLEX.out.fastq
    ch_multiqc_files    = ch_multiqc_files.mix(DEMULTIPLEX.out.reports)
    ch_versions         = ch_versions.mix(DEMULTIPLEX.out.versions)

    //
    // FASTQ QC
    //
    // MODULE: fastp
    // Run QC, trimming and adapter removal
    ch_trimmed_fastq    = FASTP(ch_raw_fastq, false, false).out.reads
    ch_multiqc_files    = ch_multiqc_files.mix( FASTP.out.json.map { meta, json -> return json} )
    ch_versions         = ch_versions.mix(FASTP.out.versions)

    //
    // ALIGNMENT
    //
    // MODULE: bowtie2
    // Align fastq files to reference genome
    ch_split_bam        = BOWTIE2(ch_trimmed_fastq, [], false).out.bam
    ch_versions         = ch_versions.mix(BOWTIE2.out.versions)

    //
    // POST ALIGNMENT
    //

    // TODO: Gather bam records per sample into ch_split_bam_per_sample

    // MODULE: biobambam/bamsormadup
    // Take multiple split bam files per sample, merge, sort and mark duplicates
    ch_merged_bam       = BAMSORMADUP(ch_split_bam_per_sample, [], false).out.bam
    ch_merged_bai       = BAMSORMADUP(ch_split_bam_per_sample, [], false).out.bam_index
    ch_merged_bam_bai   = ch_merged_bam.join(ch_merged_bai)
    ch_multiqc_files    = ch_multiqc_files.mix( BAMSORMADUP.out.metrics.map { meta, metrics -> return metrics} )
    ch_versions         = ch_versions.mix(BAMSORMADUP.out.versions)

    // MODULE: stadeniolib/scramble
    // Compress bam to cram
    ch_cram             = SCRAMBLE(ch_merged_bam, [], false).out.cram
    ch_versions         = ch_versions.mix(SCRAMBLE.out.versions)

    //
    // COVERAGE
    //
    // MODULE: mosdepth
    // Generate coverage beds
    MOSDEPTH(ch_merged_bam, [], false)
    ch_multiqc_files    = ch_multiqc_files.mix(
        MOSDEPTH.out.summary_txt.map { meta, summary_txt -> return summary_txt},
        MOSDEPTH.out.global_txt.map  { meta, global_txt -> return global_txt},
        MOSDEPTH.out.regions_txt.map { meta, regions_txt -> return regions_txt}
    )
    ch_versions         = ch_versions.mix(MOSDEPTH.out.versions)

    //
    // QC
    // SUBWORKFLOW: bam_stats_samtools
    // Run samtools QC modules
    BAM_STATS_SAMTOOLS(ch_merged_bam, [], false)
    ch_multiqc_files    = ch_multiqc_files.mix(
        BAM_STATS_SAMTOOLS.out.stats.map    { meta, stats -> return stats},
        BAM_STATS_SAMTOOLS.out.flagstat.map { meta, flagstat -> return flagstat},
        BAM_STATS_SAMTOOLS.out.idxstats.map { meta, idxstats -> return idxstats}
    )
    ch_versions         = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)

    // SUBWORKFLOW: bam_stats_picard
    // Run Picard QC modules
    BAM_QC_PICARD(ch_merged_bam, [], false)
    ch_multiqc_files    = ch_multiqc_files.mix(
        BAM_QC_PICARD.out.coverage_metrics.map { meta, coverage_metrics -> return coverage_metrics},
        BAM_QC_PICARD.out.multiple_metrics.map { meta, multiple_metrics -> return multiple_metrics},
    )
    ch_versions         = ch_versions.mix(BAM_QC_PICARD.out.versions)

    // MODULE: CUSTOM_DUMPSOFTWAREVERSIONS
    // Gather software versions for QC report
    CUSTOM_DUMPSOFTWAREVERSIONS (ch_versions.unique().collectFile(name: 'collated_versions.yml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    // MODULE: MultiQC
    // Generate aggregate QC report
    workflow_summary    = WorkflowCmggpreprocessing.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = ch_multiqc_files.mix(
        Channel.from(ch_multiqc_config),
        ch_multiqc_custom_config.collect().ifEmpty([]),
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
    )

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
