/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowCmgg-preprocessing.initialise(params, log)

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

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

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
include { BAM_TO_FASTQ} from "../subworkflows/local/bam_to_fastq/main"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { BIOBAMBAM_BAMSORMADUP       } from "../modules/nf-core/modules/biobambam/bamsormadup/main"
include { CAT_FASTQ                   } from '../modules/nf-core/modules/cat/fastq/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from "../modules/nf-core/modules/custom/dumpsoftwareversions/main"
include { FASTP                       } from "../modules/nf-core/modules/fastp/main"
include { FGBIO_FASTQTOBAM            } from "../modules/nf-core/modules/fgbio/fastqtobam/main"
include { MULTIQC                     } from "../modules/nf-core/modules/multiqc/main"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow CMGG-PREPROCESSING {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    // Sanitize inputs and separate input types
    // For now assume 1 bam/cram per sample
    ch_inputs = extract_csv(ch_input).branch {
        fastq   : (it.size() == 2 && it[1] instanceof List)
        flowcell: (it.size() == 4 && it[1].toString().endsWith(".csv"))
        reads   : (it.size() == 2 && (it[1].toString().endsWith(".bam") || it[1].toString().endsWith(".cram")))
        other   : true
    }

    ch_inputs.fastq.dump(tag: "fastq inputs", {FormattingService.prettyFormat(it)})
    ch_inputs.flowcell.dump(tag: "flowcell inputs", {FormattingService.prettyFormat(it)})
    ch_inputs.reads.dump(tag: "reads inputs", {FormattingService.prettyFormat(it)})
    ch_inputs.other.dump(tag: "other inputs", {FormattingService.prettyFormat(it)})
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // PROCESS FLOWCELL INPUTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_flowcell = ch_inputs.flowcell.multiMap { meta, samplesheet, flowcell, sample_info ->
        fc   : [meta, samplesheet, flowcell]
        info : sample_info
    }

    // DEMULTIPLEX([meta, samplesheet, flowcell])
    DEMULTIPLEX(ch_flowcell.fc)
    DEMULTIPLEX.out.bclconvert_fastq.dump(tag: "bclconvert_fastq",{FormattingService.prettyFormat(it)})
    ch_multiqc_files = ch_multiqc_files.mix(DEMULTIPLEX.out.bclconvert_reports.map { meta, reports -> return reports} )
    ch_versions      = ch_versions.mix(DEMULTIPLEX.out.versions)

    // Add metadata to demultiplexed fastq's
    ch_demultiplexed_fastq = merge_sample_info(
        DEMULTIPLEX.out.bclconvert_fastq,
        parse_sample_info_csv(ch_flowcell.info)
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // PROCESS BAM/CRAM INPUTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Convert bam/cram inputs to fastq
    BAM_TO_FASTQ(ch_inputs.reads, params.fasta)
    ch_converted_fastq = BAM_TO_FASTQ.out.fastq
    ch_converted_fastq.dump(tag: "converted_fastq",{FormattingService.prettyFormat(it)})


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // GATHER PROCESSED INPUTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // "Gather" fastq's from demultiplex and fastq inputs
    ch_sample_fastqs = Channel.empty()
    ch_sample_fastqs = ch_sample_fastqs.mix(ch_inputs.fastq, ch_demultiplexed_fastq, ch_converted_fastq)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // FASTQ TRIMMING AND QC
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // MODULE: fastp
    // Run QC, trimming and adapter removal
    // FASTP([meta, fastq], save_trimmed, save_merged)
    FASTP(ch_sample_fastqs, false, false)
    FASTP.out.reads.dump(tag: "fastp_trimmed_reads",{FormattingService.prettyFormat(it)})
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.map { meta, json -> return json} )
    ch_versions      = ch_versions.mix(FASTP.out.versions)

    // merge FASTP.out.reads in channel per sample and split samples into human and non human data
    ch_trimmed_reads = gather_split_files_per_sample(FASTP.out.reads)
        .branch { meta, reads ->
            human: meta.organism ==~ /(?i)Homo sapiens/
            other: true
        }
    ch_trimmed_reads.human.dump(tag: "human_reads",{FormattingService.prettyFormat(it)})
    ch_trimmed_reads.other.dump(tag: "other_reads",{FormattingService.prettyFormat(it)})

    // if aliger == snapaligner, proceed
    // else split into channel per chunk

    ch_reads_to_map = params.aligner == "snapaligner" ? ch_trimmed_reads.human : ch_trimmed_reads.human.map{meta, reads ->
        reads_files = meta.single_end ? reads : reads.sort().collate(2)
        return [meta, reads_files]
    }.transpose()
    ch_reads_to_map.dump(tag: "reads_to_map",{FormattingService.prettyFormat(it)})

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // STEP: FASTQ TO BAM CONVERSION
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Convert non-standard fastq data (e.g. non-human, non-DNA, ...) to BAM
    // CAT_FASTQ([meta, fastq])
    // Merge split fastqs to simplify converting to uBAM
    CAT_FASTQ(ch_trimmed_reads.other)
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)

    // FGBIO_FASTQTOBAM([meta, fastq])
    FGBIO_FASTQTOBAM(CAT_FASTQ.out.reads)
    ch_versions = ch_versions.mix(FGBIO_FASTQTOBAM.out.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // STEP: ALIGNMENT
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // align fastq files per sample, merge, sort and markdup.
    // ALIGNMENT([meta,fastq], index, sort)
    ALIGNMENT(ch_reads_to_map, ch_map_index, true)
    ch_versions = ch_versions.mix(ALIGNMENT.out.versions)

    // Gather bams per sample for merging
    ch_bam_per_sample = params.aligner == "snapaligner" ? ALIGNMENT.out.bam : gather_split_files_per_sample(ALIGNMENT.out.bam)
    ch_bam_per_sample.dump(tag: "bam_per_sample",{FormattingService.prettyFormat(it)})

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // STEP: MARK DUPLICATES
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // BIOBAMBAM_BAMSORMADUP([meta, [bam, bam]], fasta)
    // Should only run when alinger != snapaligner
    BIOBAMBAM_BAMSORMADUP(ch_bam_per_sample, params.fasta)
    ch_markdup_bam_bai = BIOBAMBAM_BAMSORMADUP.out.bam.join(BIOBAMBAM_BAMSORMADUP.out.bam_index)
    ch_multiqc_files = ch_multiqc_files.mix( BIOBAMBAM_BAMSORMADUP.out.metrics.map { meta, metrics -> return metrics} )
    ch_versions = ch_versions.mix(BIOBAMBAM_BAMSORMADUP.out.versions)

    ch_markdup_bam_bai = Channel.empty()
    ch_markdup_bam_bai = ch_markdup_bam_bai.mix(ALIGNMENT.out.bam.join(ALIGNMENT.out.bai), BIOBAMBAM_BAMSORMADUP.out.bam.join(BIOBAMBAM_BAMSORMADUP.out.bam_index))
    ch_markdup_bam_bai.dump(tag: "markdup_bam_bai",{FormattingService.prettyFormat(it)})

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // COVERAGE ANALYSIS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Generate coverage metrics and beds for each sample
    // COVERAGE([meta,bam, bai], fasta, fai, target, bait)
    COVERAGE(
        ch_markdup_bam_bai,
        params.fasta, params.fai,
        [],
        []
    )
    ch_multiqc_files = ch_multiqc_files.mix( COVERAGE.out.metrics.map { meta, metrics -> return metrics} )
    ch_versions      = ch_versions.mix(COVERAGE.out.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // QC FOR BAM FILES
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Gather metrics from bam files
    BAM_QC(
        ch_markdup_bam_bai,   // [meta, bam, bai]
    )
    ch_multiqc_files = ch_multiqc_files.mix( BAM_QC.out.metrics.map { meta, metrics -> return metrics} )
    ch_versions      = ch_versions.mix(BAM_QC.out.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // ARCHIVING
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Compress and checksum bam files
    ch_reads_to_compress = Channel.empty()
    ch_reads_to_compress = ch_reads_to_compress.mix(
        ch_markdup_bam_bai,
        FGBIO_FASTQTOBAM.out.bam.map({ meta, bam ->
            new_meta = meta.clone()
            new_meta.id = "${meta.id}.unaligned"
            return [ new_meta, bam, [] ]
        })
    )
    ch_reads_to_compress.dump(tag: "reads_to_compress", {FormattingService.prettyFormat(it)})

    BAM_ARCHIVE(
        ch_reads_to_compress,
        params.fasta,
        params.fai
    )
    ch_versions = ch_versions.mix(BAM_ARCHIVE.out.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // REPORTING
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // MODULE: CUSTOM_DUMPSOFTWAREVERSIONS
    // Gather software versions for QC report
    CUSTOM_DUMPSOFTWAREVERSIONS (ch_versions.unique().collectFile(name: "collated_versions.yml"))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    // MODULE: MULTIQC
    // Generate aggregate QC report
    workflow_summary    = WorkflowCmggpreprocessing.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowCmgg-preprocessing.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.collect().ifEmpty([]),
        ch_multiqc_custom_config.collect().ifEmpty([]),
        ch_multiqc_logo.collect().ifEmpty([])
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Extract information (meta data + file(s)) from csv file(s)
def extract_csv(csv_file) {
    Channel.value(csv_file).splitCsv(header: true, strip: true).map { row ->
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
            log.error "Fastq input requires both samplename and fastq"
        }
        // check for invalid bam/cram input
        if((row.bam || row.cram) && !(row.samplename)){
            log.error "BAM/CRAM input requires both samplename and bam/cram"
        }
        // check for mixed input
        if(row.flowcell && row.fastq_1){
            log.error "Cannot mix flowcell and fastq/cram inputs"
        }
        // valid flowcell input
        if(row.flowcell && row.samplesheet){
            return parse_flowcell_csv(row)
        }
        // valid fastq input
        if(row.fastq_1 && row.samplename){
            return parse_fastq_csv(row)
        }
        // valid bam/cram input
        if((row.bam || row.cram) && row.samplename){
            return parse_reads_csv(row)
        }
    }
}

// Parse flowcell input map
def parse_flowcell_csv(row) {
    def meta = [:]
    meta.id   = row.id.toString()
    meta.lane = null
    if (row.containsKey("lane") && row.lane ) {
        meta.lane = row.lane.toInteger()
    }

    def flowcell        = file(row.flowcell, checkIfExists: true)
    def samplesheet     = file(row.samplesheet, checkIfExists: true)
    if (!(row.sample_info)) log.error "Sample info csv not defined"
    def sample_info_csv = file(row.sample_info, checkIfExists: true)
    return [meta, samplesheet, flowcell, sample_info_csv]
}

// Parse fastq input map
def parse_fastq_csv(row) {
    def fastq_1 = file(row.fastq_1, checkIfExists: true)
    def fastq_2 = row.fastq_2 ? file(row.fastq_2, checkIfExists: true) : null

    def meta = [:]
    meta.id         = row.id.toString()
    meta.samplename = row.samplename.toString()
    meta.organism   = row.organism ? row.organism.toString() : ""
    meta.single_end = fastq_2      ? false : true
    // Set readgroup info
    meta.readgroup    = [:]
    meta.readgroup    = readgroup_from_fastq(fastq_1)
    meta.readgroup.SM = meta.samplename
    meta.readgroup.LB = row.library ? row.library.toString() : ""
    meta.readgroup.ID = meta.readgroup.ID ? meta.readgroup.ID : meta.samplename

    return [meta, fastq_2 ? [fastq_1, fastq_2] : [fastq_1]]
}

// Parse sample info input map
def parse_sample_info_csv(csv_file) {
    csv_file.splitCsv(header: true, strip: true).map { row ->
        // check mandatory fields
        if (!(row.samplename)) log.error "Missing samplename field in sample info file"
        return row
    }
}

// Parse bam/cram input map
def parse_reads_csv(row) {
    def cram = row.cram ? file(row.cram, checkIfExists: true) : null

    def bam = row.bam ? file(row.bam, checkIfExists: true) : null

    def meta = [:]
    meta.id         = row.id.toString()
    meta.samplename = row.samplename.toString()
    meta.organism   = row.organism ? row.organism.toString() : ""

    return [meta, bam ? bam : cram]
}

// Merge fastq meta with sample info
def merge_sample_info(ch_fastq, ch_sample_info) {
    ch_fastq
    .combine(ch_sample_info)
    .map { meta1, fastq, meta2 ->
        def meta = meta1.clone()
        if ( meta2 && (meta1.samplename == meta2.samplename)) {
            meta = meta1 + meta2
            meta.readgroup    = [:]
            meta.readgroup    = readgroup_from_fastq(fastq[0])
            meta.readgroup.SM = meta.samplename
            if(meta.library){
                meta.readgroup.LB =  meta.library.toString()
            }
            return [ meta, fastq ]
        }
    }.groupTuple( by: [0])
    .map { meta, fq ->
        return [meta, fq.flatten().unique()]
    }
}

// https://github.com/nf-core/sarek/blob/7ba61bde8e4f3b1932118993c766ed33b5da465e/workflows/sarek.nf#L1014-L1040
def readgroup_from_fastq(path) {
    // expected format:
    // xx:yy:FLOWCELLID:LANE:... (seven fields)
    // or
    // FLOWCELLID:LANE:xx:... (five fields)
    def line

    path.withInputStream {
        InputStream gzipStream = new java.util.zip.GZIPInputStream(it)
        Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
        BufferedReader buffered = new BufferedReader(decoder)
        line = buffered.readLine()
    }
    assert line.startsWith('@')
    line = line.substring(1)
    def fields = line.split(':')
    def rg = [:]
    rg.CN = "CMGG"

    if (fields.size() >= 7) {
        // CASAVA 1.8+ format, from  https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.htm
        // "@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>:<UMI> <read>:<is filtered>:<control number>:<index>"
        sequencer_serial = fields[0]
        run_nubmer       = fields[1]
        fcid             = fields[2]
        lane             = fields[3]
        index            = fields[-1] =~ /[GATC+-]/ ? fields[-1] : ""

        rg.ID = fcid
        rg.PU = [fcid, index].findAll().join(".")
        rg.PL = "ILLUMINA"
    } else if (fields.size() == 5) {
        fcid = fields[0]
        rg.ID = fcid
    }
    return rg
}

// Gather split files per sample
def gather_split_files_per_sample(ch_files) {
    // Gather bam files per sample based on id
    ch_files.map {
        // set id to filename without lane designation
        meta, files ->
        new_meta = meta.clone()
        new_meta.id = meta.samplename
        return [new_meta, files]
    }
    .groupTuple( by: [0])
    .map { meta, files ->
        return [meta, files.flatten()]
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
    if (params.hook_url) {
        NfcoreTemplate.adaptivecard(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.adaptivecard(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
