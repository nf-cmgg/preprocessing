/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowCmggpreprocessing.initialise(params, log)

def checkPathParamList = [ params.input, params.samples, params.multiqc_config, params.fasta, params.fasta_fai ]
// TODO: restore this
//for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, "Flowcell metadata sheet not specified!" }
if (params.samples) { ch_input = file(params.samples) } else { exit 1, "Sample metadata sheet not specified!" }


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
include { BAM_STATS_SAMTOOLS} from "../subworkflows/nf-core/subworkflows/bam_stats_samtools/main"
include { DEMULTIPLEX       } from "../subworkflows/local/demultiplex"
include { INPUT_CHECK       } from "../subworkflows/local/input_check"
//inlcude { BAM_QC_PICARD     } from "../subworkflows/nf-core/subworkflows/bam_qc_picard/main"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { BIOBAMBAM_BAMSORMADUP       } from "../modules/nf-core/modules/biobambam/bamsormadup/main"
include { BOWTIE2_ALIGN               } from "../modules/nf-core/modules/bowtie2/align/main"
include { CUSTOM_DUMPSOFTWAREVERSIONS } from "../modules/nf-core/modules/custom/dumpsoftwareversions/main"
include { FASTP                       } from "../modules/nf-core/modules/fastp/main"
include { MD5SUM                      } from "../modules/nf-core/modules/md5sum/main"
include { MOSDEPTH                    } from "../modules/nf-core/modules/mosdepth/main"
include { MULTIQC                     } from "../modules/nf-core/modules/multiqc/main"
include { SAMTOOLS_INDEX              } from "../modules/nf-core/modules/samtools/index/main"
include { SAMTOOLS_CONVERT            } from "../modules/nf-core/modules/samtools/convert/main"

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

    //*
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //*
    ch_flowcells        = INPUT_CHECK (ch_input).flowcells
    ch_versions         = ch_versions.mix(INPUT_CHECK.out.versions)

    //*
    // DEMULTIPLEXING
    //*
    // SUBWORKFLOW: demultiplex
    DEMULTIPLEX(ch_flowcells)
    ch_multiqc_files    = ch_multiqc_files.mix(DEMULTIPLEX.out.reports)
    ch_versions         = ch_versions.mix(DEMULTIPLEX.out.versions)

    //*
    // FASTQ QC and TRIMMING
    //*

    // MODULE: fastp
    // Run QC, trimming and adapter removal
    FASTP(DEMULTIPLEX.out.fastq, false, false)
    ch_multiqc_files    = ch_multiqc_files.mix( FASTP.out.json.map { meta, json -> return json} )
    ch_versions         = ch_versions.mix(FASTP.out.versions)

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
