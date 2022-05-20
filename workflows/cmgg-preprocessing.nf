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
include { BAM_QC_PICARD         } from "../subworkflows/nf-core/subworkflows/bam_qc_picard/main"
include { BAM_STATS_SAMTOOLS    } from "../subworkflows/nf-core/subworkflows/bam_stats_samtools/main"
include { MARKDUP_PARALLEL      } from "../subworkflows/local/markdup_parallel/main"
include { DEMULTIPLEX           } from "../subworkflows/local/demultiplex/main"
include { INPUT_CHECK           } from "../subworkflows/local/input_check"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { BOWTIE2_ALIGN               } from "../modules/nf-core/modules/bowtie2/align/main"
include { CUSTOM_DUMPSOFTWAREVERSIONS } from "../modules/nf-core/modules/custom/dumpsoftwareversions/main"
include { FASTP                       } from "../modules/nf-core/modules/fastp/main"
include { MD5SUM                      } from "../modules/nf-core/modules/md5sum/main"
include { MOSDEPTH                    } from "../modules/nf-core/modules/mosdepth/main"
include { MULTIQC                     } from "../modules/nf-core/modules/multiqc/main"
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
    // STEP: DEMULTIPLEXING
    //*
    // SUBWORKFLOW: demultiplex
    DEMULTIPLEX(ch_flowcells)
    ch_multiqc_files    = ch_multiqc_files.mix(DEMULTIPLEX.out.reports)
    ch_versions         = ch_versions.mix(DEMULTIPLEX.out.versions)

    //*
    // STEP: FASTQ QC and TRIMMING
    //*

    // MODULE: fastp
    // Run QC, trimming and adapter removal
    FASTP(DEMULTIPLEX.out.fastq, false, false)
    ch_multiqc_files    = ch_multiqc_files.mix( FASTP.out.json.map { meta, json -> return json} )
    ch_versions         = ch_versions.mix(FASTP.out.versions)

    //*
    // STEP: ALIGNMENT
    //*

    // MODULE: bowtie2
    // Align fastq files to reference genome and sort
    // BOWTIE2_ALIGN([meta,reads, index, save_unaligned, sort])
    BOWTIE2_ALIGN(FASTP.out.reads, params.bowtie2, false, true)
    ch_versions         = ch_versions.mix(BOWTIE2_ALIGN.out.versions)

    //*
    // STEP: POST ALIGNMENT
    //*

    // gather split bams per samples
    ch_bam_per_sample = gather_bam_per_sample(BOWTIE2_ALIGN.out.bam)

    // SUBWORKFLOW: parallel markdup
    // Take split bam files, mark duplicates and merge
    MARKDUP_PARALLEL(ch_bam_per_sample)


    // MODULE: samtools/convert
    // Compress bam to cram
    SAMTOOLS_CONVERT(
        MARKDUP_PARALLEL.out.bam_bai.map {
            meta, bam, bai -> return [meta, bam]
        }, params.fasta, params.fai
    )
    ch_versions         = ch_versions.mix(SAMTOOLS_CONVERT.out.versions)

    // MODULE: MD5SUM
    // Generate md5sum for cram file
    MD5SUM(
        SAMTOOLS_CONVERT.out.alignment_index.map {
            meta, cram, crai -> return [meta, cram]
        }
    )
    ch_versions         = ch_versions.mix(MD5SUM.out.versions)

    //*
    // COVERAGE
    //*

    // MODULE: mosdepth
    // Generate coverage beds
    MOSDEPTH(MARKDUP_PARALLEL.out.bam_bai, [], false)
    ch_multiqc_files    = ch_multiqc_files.mix(
        MOSDEPTH.out.summary_txt.map { meta, summary_txt -> return summary_txt},
        MOSDEPTH.out.global_txt.map  { meta, global_txt -> return global_txt},
        MOSDEPTH.out.regions_txt.map { meta, regions_txt -> return regions_txt}
    )
    ch_versions         = ch_versions.mix(MOSDEPTH.out.versions)

    //*
    // QC
    //*

    // SUBWORKFLOW: bam_stats_samtools
    // Run samtools QC modules
    BAM_STATS_SAMTOOLS(MARKDUP_PARALLEL.out.bam_bai)
    ch_multiqc_files    = ch_multiqc_files.mix(
        BAM_STATS_SAMTOOLS.out.stats.map    { meta, stats -> return stats},
        BAM_STATS_SAMTOOLS.out.flagstat.map { meta, flagstat -> return flagstat},
        BAM_STATS_SAMTOOLS.out.idxstats.map { meta, idxstats -> return idxstats}
    )
    ch_versions         = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)

    // SUBWORKFLOW: bam_stats_picard
    // Run Picard QC modules
    BAM_QC_PICARD(
        MARKDUP_PARALLEL.out.bam_bai.map {
            meta, bam, bai -> return [meta, bam]
        },
        [],
        false
    )
    ch_multiqc_files    = ch_multiqc_files.mix(
        BAM_QC_PICARD.out.coverage_metrics.map { meta, coverage_metrics -> return coverage_metrics},
        BAM_QC_PICARD.out.multiple_metrics.map { meta, multiple_metrics -> return multiple_metrics},
    )
    ch_versions         = ch_versions.mix(BAM_QC_PICARD.out.versions)


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
