/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowCmggpreprocessing.initialise(params, log)

def checkPathParamList = [ params.input, params.samples, params.multiqc_config, params.fasta, params.fai ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, "Flowcell metadata sheet not specified!" }
if (params.samples) { ch_input = file(params.samples) } else { exit 1, "Sample metadata sheet not specified!" }


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

multiqc_config = params.multiqc_config ? file(params.multiqc_config, checkIfExists: true) : file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
multiqc_logo   = params.multiqc_logo ? file(params.multiqc_logo, checkIfExists: true) : file("$projectDir/assets/CMGG_logo.png", checkIfExists: true)

// Info required for completion email and summary
def multiqc_report = []

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from "../subworkflows/local/input_check"
include { DEMULTIPLEX } from "../subworkflows/local/demultiplex/main"
include { ALIGNMENT   } from "../subworkflows/local/alignment/main"
include { COVERAGE    } from "../subworkflows/local/coverage/main"
include { BAM_QC      } from "../subworkflows/local/bam_qc/main"
include { BAM_ARCHIVE } from "../subworkflows/local/bam_archive/main"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS } from "../modules/nf-core/modules/custom/dumpsoftwareversions/main"
include { MULTIQC                     } from "../modules/nf-core/modules/multiqc/main"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CMGGPREPROCESSING {

    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()


    // Gather index for mapping given the chosen aligner
    ch_map_index =  params.aligner == "bwa"         ? params.bwa         :
                    params.aligner == "bwamem2"     ? params.bwamem2     :
                    params.aligner == "bowtie2"     ? params.bowtie2     :
                    params.aligner == "dragmap"     ? params.dragmap     :
                    params.aligner == "snapaligner" ? params.snapaligner :
                    []

    // Sanitize inputs
    ch_flowcells = INPUT_CHECK (ch_input).flowcells
    ch_versions  = ch_versions.mix(INPUT_CHECK.out.versions)

    //*
    // STEP: DEMULTIPLEXING and FASTQ QC
    //*
    // DEMULTIPLEX([meta, samplesheet, flowcell])
    DEMULTIPLEX(ch_flowcells)
    ch_multiqc_files = ch_multiqc_files.mix(
        DEMULTIPLEX.out.bclconvert_reports.map { meta, reports -> return reports},
        DEMULTIPLEX.out.fastp_reports.map { meta, json -> return json}
    )
    ch_versions = ch_versions.mix(DEMULTIPLEX.out.versions)

    // TODO: parse metadata from params.samples and merge with metadata from DEMULTIPLEX.out.fastq
    // TODO: Unblock workflow by first collecting ALL fastq's per sample, then running subworkflow per sample
    // Currently, the wf is blocked by the alignment step, which makes all other steps wait for the alignment step to finish
    // "gather_bam_per_sample" is the culprit here

    ch_fastq_per_sample = DEMULTIPLEX.out.bclconvert_fastq

    //*
    // STEP: ALIGNMENT
    //*
    // align fastq files per sample, merge, sort and markdup.
    // ALIGNMENT([meta,fastq], index, sort)
    ALIGNMENT(ch_fastq_per_sample, ch_map_index, true)
    ch_multiqc_files = ch_multiqc_files.mix(ALIGNMENT.out.markdup_metrics.map { meta, metrics -> return metrics})
    ch_versions = ch_versions.mix(ALIGNMENT.out.versions)

    //*
    // STEP: COVERAGE
    //*
    // Generate coverage metrics and beds for each sample
    // COVERAGE([meta,bam, bai], fasta, fai, target, bait)
    COVERAGE(
        ALIGNMENT.out.bam_bai,
        params.fasta, params.fai,
        [],
        []
    )
    ch_multiqc_files    = ch_multiqc_files.mix( COVERAGE.out.metrics.map { meta, metrics -> return metrics} )
    ch_versions = ch_versions.mix(COVERAGE.out.versions)

    //*
    // STEP: QC
    //*
    // Gather metrics from bam files
    BAM_QC(
        ALIGNMENT.out.bam_bai,   // [meta, bam, bai]
    )
    ch_multiqc_files    = ch_multiqc_files.mix( BAM_QC.out.metrics.map { meta, metrics -> return metrics} )
    ch_versions = ch_versions.mix(BAM_QC.out.versions)

    //*
    // STEP: BAM ARCHIVE
    //*
    // Compress and checksum bam files
    BAM_ARCHIVE(
        ALIGNMENT.out.bam_bai.map {meta, bam, bai -> return [meta,bam]},
        params.fasta,
        params.fai
    )

    //*
    // STEP: POST QC
    //*
    // MODULE: CUSTOM_DUMPSOFTWAREVERSIONS
    // Gather software versions for QC report
    CUSTOM_DUMPSOFTWAREVERSIONS (ch_versions.unique().collectFile(name: "collated_versions.yml"))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    // MODULE: MULTIQC
    // Generate aggregate QC report
    workflow_summary    = WorkflowCmggpreprocessing.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml"),
    )

    MULTIQC (
        ch_multiqc_files.collect(), [multiqc_config, multiqc_logo]
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def gather_bam_per_sample(ch_aligned_bam) {
    // Gather bam files per sample based on id
    ch_aligned_bam.map {
        // set id to filename without lane designation
        meta, bam ->
        new_meta = meta.clone()
        new_meta.id = meta.id - ~/_S[0-9]+_.*$/
        return [new_meta, bam]
    }.groupTuple( by: [0])
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
