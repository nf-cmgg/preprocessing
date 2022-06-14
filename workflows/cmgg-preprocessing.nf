/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowCmggpreprocessing.initialise(params, log)

def checkPathParamList = [ params.input, params.multiqc_config, params.fasta, params.fai ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, "inputs sheet not specified!" }


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
include { BIOBAMBAM_BAMSORMADUP       } from "../modules/nf-core/modules/biobambam/bamsormadup/main"
include { CUSTOM_DUMPSOFTWAREVERSIONS } from "../modules/nf-core/modules/custom/dumpsoftwareversions/main"
include { FASTP                       } from "../modules/nf-core/modules/fastp/main"
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

    // Sanitize inputs and separate input types
    ch_inputs = extract_csv(ch_input).branch {
        fastq: it.size() == 2
        flowcell: it.size() == 4
    }

    //*
    // STEP: DEMULTIPLEX FLOWCELLS
    //*
    // DEMULTIPLEX([meta, samplesheet, flowcell])
    DEMULTIPLEX(
        ch_inputs.flowcell.map {
            meta,samplesheet, flowcell, sample_info ->
            return [meta, samplesheet, flowcell]
        }
    )
    ch_multiqc_files = ch_multiqc_files.mix(
        DEMULTIPLEX.out.bclconvert_reports.map { meta, reports -> return reports},
    )
    ch_versions = ch_versions.mix(DEMULTIPLEX.out.versions)

    //*
    // STEP: FASTQ TRIMMING AND QC
    //*
    // "Gather" fastq's from demultiplex and fastq inputs
    ch_sample_fastqs = Channel.empty()
    ch_sample_fastqs = ch_sample_fastqs.mix(
        ch_inputs.fastq
        DEMULTIPLEX.out.bclconvert_fastq,
    )

    // MODULE: fastp
    // Run QC, trimming and adapter removal
    // FASTP([meta, fastq], save_trimmed, save_merged)
    FASTP(ch_sample_fastqs, false, false)
    ch_multiqc_files  = ch_multiqc_files.mix(
        FASTP.out.json.map { meta, json -> return json}
    )
    ch_versions = ch_versions.mix(FASTP.out.versions)

    //*
    // STEP: ALIGNMENT
    //*
    // align fastq files per sample, merge, sort and markdup.
    // ALIGNMENT([meta,fastq], index, sort)
    ALIGNMENT(FASTP.out.reads, ch_map_index, true)
    ch_versions = ch_versions.mix(ALIGNMENT.out.versions)

    // Gather bams per sample for merging
    ch_bam_per_sample = gather_bam_per_sample(ALIGNMENT.out.bam)

    //*
    // STEP: MARK DUPLICATES
    //*
    // BIOBAMBAM_BAMSORMADUP([meta, [bam, bam]], fasta)
    BIOBAMBAM_BAMSORMADUP(ch_bam_per_sample, params.fasta)
    ch_markdup_bam_bai = BIOBAMBAM_BAMSORMADUP.out.bam.join(BIOBAMBAM_BAMSORMADUP.out.bam_index)
    ch_multiqc_files = ch_multiqc_files.mix( BIOBAMBAM_BAMSORMADUP.out.metrics.map { meta, metrics -> return metrics} )
    ch_versions = ch_versions.mix(BIOBAMBAM_BAMSORMADUP.out.versions)

    //*
    // STEP: COVERAGE
    //*
    // Generate coverage metrics and beds for each sample
    // COVERAGE([meta,bam, bai], fasta, fai, target, bait)
    COVERAGE(
        ch_markdup_bam_bai,
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
        ch_markdup_bam_bai,   // [meta, bam, bai]
    )
    ch_multiqc_files    = ch_multiqc_files.mix( BAM_QC.out.metrics.map { meta, metrics -> return metrics} )
    ch_versions = ch_versions.mix(BAM_QC.out.versions)

    //*
    // STEP: BAM ARCHIVE
    //*
    // Compress and checksum bam files
    BAM_ARCHIVE(
        ch_markdup_bam_bai.map {meta, bam, bai -> return [meta,bam]},
        params.fasta,
        params.fai
    )
    ch_versions = ch_versions.mix(BAM_ARCHIVE.out.versions)

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

// Function to extract information (meta data + file(s)) from csv file(s)
def extract_csv(csv_file) {
    Channel.from(csv_file).splitCsv(header: true, strip: true).map { row ->
        // check common mandatory fields
        if(!(row.id)){
            log.error "Missing id field in input csv file"
        }
        // check for invalid flowcell input
        if(row.flowcell && !(row.samplesheet)){
            log.error "Flowcell input requires both samplesheet and flowcell"
        }
        // check for invalid fastq input
        if(row.fastq_1 && !(row.samplename)){
            log.error "Flowcell input requires both samplesheet and flowcell"
        }
        // check for mixed input
        if(row.flowcell && row.fastq_1){
            log.error "Cannot mix flowcell and fastq inputs"
        }
        // valid flowcell input
        if(row.flowcell && row.samplesheet){
            return parse_flowcell_csv(row)
        }
        // valid fastq input
        if(row.fastq_1 && row.samplename){
            return parse_fastq_csv(row)
        }
    }
}

def parse_flowcell_csv(row) {
    def meta = [:]

    meta.id = row.id
    meta.lane = row.lane ?: ""

    def flowcell    = file(row.flowcell, checkIfExists: true)
    def samplesheet = file(row.samplesheet, checkIfExists: true)
    def sample_info = row.sample_info ? file(row.sample_info, checkIfExists: true) : null

    return [meta, samplesheet, flowcell, sample_info]
}

def parse_fastq_csv(row) {
    def meta = [:]

    meta.id         = row.id
    meta.samplename = row.samplename
    meta.readgroup  = row.readgroup ? row.readgroup.toString() : ""
    meta.reference  = row.reference ? row.reference.toString() : ""
    meta.single_end = row.fastq_2   ? false : true

    def fastq_1 = file(row.fastq_1, checkIfExists: true)
    def fastq_2 = row.fastq_2 ? file(row.fastq_2, checkIfExists: true) : null
    def reads = fastq_2 ? [fastq_1, fastq_2] : [fastq_1]
    return [meta, reads]
}

def parse_sample_info_csv(csv_file) {
    Channel.from(csv_file).splitCsv(header: true).map { row ->
        // check mandatory fields
        if (!(row.id && row.samplename)){
            log.error "Missing id or samplename field in sample info file"
        } else {
            return [row]
        }
    }
}

// Function to gather bam files per sample
def gather_bam_per_sample(ch_aligned_bam) {
    // Gather bam files per sample based on id
    ch_aligned_bam.map {
        // set id to filename without lane designation
        meta, bam ->
        new_meta = meta.clone()
        new_meta.id = meta.id - ~/_S[0-9]+_.*$/
        return [new_meta, bam]
    }
    .groupTuple( by: [0])
    .map { meta, bam ->
        return [meta, bam.flatten()]
    }
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
