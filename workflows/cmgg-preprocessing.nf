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
include { BCL_DEMULTIPLEX   } from "../subworkflows/nf-core/bcl_demultiplex/main"
include { FASTA_INDEX_DNA   } from "../subworkflows/nf-core/fasta_index_dna/main"
include { FASTQ_TO_CRAM     } from "../subworkflows/local/fastq_to_cram/main"
include { FASTQ_TO_UCRAM    } from "../subworkflows/local/fastq_to_ucram/main"
include { COVERAGE          } from "../subworkflows/local/coverage/main"
include { BAM_QC            } from "../subworkflows/local/bam_qc/main"
include { BAM_TO_FASTQ      } from "../subworkflows/local/bam_to_fastq/main"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS } from "../modules/nf-core/custom/dumpsoftwareversions/main"
include { FASTP                       } from "../modules/nf-core/fastp/main"
include { MULTIQC                     } from "../modules/nf-core/multiqc/main"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CMGGPREPROCESSING {

    // input values
    aligner = params.aligner
    postprocessor = params.postprocessor
    genome  = params.genome

    // input options
    run_coverage   = params.run_coverage
    disable_picard = params.disable_picard_metrics

    // input channels
    ch_fasta     = Channel.value([
        [id:genome],
        file(params.fasta, checkIfExists: true)
    ])
    ch_fasta_fai = Channel.value([
        [id:genome],
        file(params.fasta, checkIfExists: true),
        file(params.fai, checkIfExists: true)
    ])
    //TODO: fix this input
    ch_altliftover = Channel.value([
        [id:genome],
        []
    ])

    ch_dict  = Channel.value([
        [id:genome],
        file(params.dict,  checkIfExists: true)
    ])

    ch_aligner_index = Channel.empty()

    ch_bait_regions   = params.bait_regions   ? Channel.value(file(params.bait_regions, checkIfExists: true))   : []
    ch_target_regions = params.target_regions ? Channel.value(file(params.target_regions, checkIfExists: true)) : []

    // output channels
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Gather index for mapping given the chosen aligner
    aligner_index = aligner == "bowtie2" ? params.bowtie2 :
                    aligner == "bwamem"  ? params.bwamem  :
                    aligner == "bwamem2" ? params.bwamem2 :
                    aligner == "dragmap" ? params.dragmap :
                    aligner == "snap"    ? params.snap    :
                    []
    if (aligner_index) {
        ch_aligner_index = Channel.value([[id:genome],file(aligner_index, checkIfExists: true)])
    } else {
        FASTA_INDEX_DNA(
            ch_fasta,
            ch_altliftover,
            aligner
        )
        ch_aligner_index = FASTA_INDEX_DNA.out.index
        ch_versions      = ch_versions.mix(FASTA_INDEX_DNA.out.versions)
        ch_aligner_index.dump(tag: "MAIN: aligner index" , {FormattingService.prettyFormat(it)})
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // INPUT PARSING
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Sanitize inputs and separate input types
    // For now assume 1 bam/cram per sample
    ch_inputs = extract_csv(ch_input).branch {
        fastq   : (it.size() == 2 && it[1] instanceof List)
        flowcell: (it.size() == 4 && it[1].toString().endsWith(".csv"))
        reads   : (it.size() == 2 && (it[1].toString().endsWith(".bam") || it[1].toString().endsWith(".cram")))
        other   : true
    }

    ch_inputs.fastq
        .dump(tag: "MAIN: fastq inputs"   , {FormattingService.prettyFormat(it)})
    ch_inputs.flowcell
        .dump(tag: "MAIN: flowcell inputs", {FormattingService.prettyFormat(it)})
    ch_inputs.reads
        .dump(tag: "MAIN: reads inputs"   , {FormattingService.prettyFormat(it)})
    ch_inputs.other
        .dump(tag: "MAIN: other inputs"   , {FormattingService.prettyFormat(it)})

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // PROCESS FLOWCELL INPUTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_flowcell = ch_inputs.flowcell.multiMap { meta, samplesheet, flowcell, sample_info ->
        fc   : [meta, samplesheet, flowcell]
        info : sample_info
    }

    // BCL_DEMULTIPLEX([meta, samplesheet, flowcell], demultiplexer)
    BCL_DEMULTIPLEX(ch_flowcell.fc, "bclconvert")
    BCL_DEMULTIPLEX.out.fastq.dump(tag: "DEMULTIPLEX: fastq",{FormattingService.prettyFormat(it)})
    ch_multiqc_files = ch_multiqc_files.mix(
        BCL_DEMULTIPLEX.out.reports.map { meta, reports -> return reports},
        BCL_DEMULTIPLEX.out.stats.map   { meta, stats   -> return stats  }
    )
    ch_versions = ch_versions.mix(BCL_DEMULTIPLEX.out.versions)

    // Add metadata to demultiplexed fastq's
    ch_demultiplexed_fastq = merge_sample_info(
        BCL_DEMULTIPLEX.out.fastq,
        parse_sample_info_csv(ch_flowcell.info)
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // PROCESS BAM/CRAM INPUTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Convert bam/cram inputs to fastq
    BAM_TO_FASTQ(ch_inputs.reads, ch_fasta_fai)
    ch_converted_fastq = BAM_TO_FASTQ.out.fastq
    ch_converted_fastq.dump(tag: "converted_fastq",{FormattingService.prettyFormat(it)})


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // GATHER PROCESSED INPUTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // "Gather" fastq's from demultiplex and fastq inputs
    ch_sample_fastqs = Channel.empty()
    ch_sample_fastqs = count_samples(
        ch_sample_fastqs.mix(
                ch_inputs.fastq, ch_demultiplexed_fastq, ch_converted_fastq
        )
    ).map { meta, reads ->
        return [meta - meta.subMap('fcid', 'lane', 'library'), reads]
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // FASTQ TRIMMING AND QC
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // MODULE: fastp
    // Run QC, trimming and adapter removal
    // FASTP([meta, fastq], save_trimmed, save_merged)
    FASTP(ch_sample_fastqs, [], false, false)
    FASTP.out.reads.dump(tag: "MAIN: fastp trimmed reads",{FormattingService.prettyFormat(it)})
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.map { meta, json -> return json} )
    ch_versions      = ch_versions.mix(FASTP.out.versions)

    // edit meta.id to match sample name
    ch_trimmed_reads = FASTP.out.reads
    .map { meta, reads ->
        read_files = meta.single_end ? reads : reads.sort{ a,b -> a.getName().tokenize('.')[0] <=> b.getName().tokenize('.')[0] }.collate(2)
        return [
            meta + [ chunks: read_files instanceof List ? read_files.size() : [read_files].size() ],
            read_files
        ]
    }
    // transpose to get read pairs
    .transpose()
    // split samples into human and non human data
    .branch { meta, reads ->
        human: meta.organism ==~ /(?i)Homo sapiens/
        other: true
    }
    ch_trimmed_reads.human.dump(tag: "MAIN: human reads",{FormattingService.prettyFormat(it)})
    ch_trimmed_reads.other.dump(tag: "MAIN: other reads",{FormattingService.prettyFormat(it)})

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // STEP: FASTQ TO UNALIGNED CRAM CONVERSION
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    FASTQ_TO_UCRAM(ch_trimmed_reads.other, ch_fasta_fai)
    ch_versions = ch_versions.mix(FASTQ_TO_UCRAM.out.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // STEP: FASTQ TO ALIGNED CRAM CONVERSION
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    FASTQ_TO_CRAM(
        ch_trimmed_reads.human.map{ meta, reads -> [meta - meta.subMap('organism'), reads]},
        ch_fasta_fai, aligner,
        ch_aligner_index, postprocessor
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_TO_CRAM.out.multiqc_files)
    ch_versions = ch_versions.mix(FASTQ_TO_CRAM.out.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // COVERAGE ANALYSIS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Generate coverage metrics and beds for each sample
    // COVERAGE([meta,bam, bai], [meta2, fasta, fai], target)
    if (run_coverage){
        COVERAGE(
            FASTQ_TO_CRAM.out.cram_crai,
            ch_fasta_fai,
            ch_target_regions,
        )
        ch_coverage_beds = Channel.empty().mix(
            COVERAGE.out.per_base_bed.join(COVERAGE.out.per_base_bed_csi),
            COVERAGE.out.regions_bed_csi.join(COVERAGE.out.regions_bed_csi),
            COVERAGE.out.quantized_bed.join(COVERAGE.out.quantized_bed_csi),
        )
        ch_multiqc_files = ch_multiqc_files.mix( COVERAGE.out.metrics.map { meta, metrics -> return metrics} )
        ch_versions      = ch_versions.mix(COVERAGE.out.versions)
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // QC FOR ALIGNMENTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Gather metrics from bam files
    // BAM_QC([meta, bam, bai], [meta2, fasta, fai], [meta2, dict], [target], [bait])
    BAM_QC( FASTQ_TO_CRAM.out.cram_crai, ch_fasta_fai, ch_dict, ch_target_regions, ch_bait_regions, disable_picard)
    ch_multiqc_files = ch_multiqc_files.mix( BAM_QC.out.metrics.map { meta, metrics -> return metrics} )
    ch_versions      = ch_versions.mix(BAM_QC.out.versions)


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

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml")
    ).collect()
    ch_multiqc_files.dump(tag: "MAIN: multiqc files",{FormattingService.prettyFormat(it)})

    MULTIQC (
        ch_multiqc_files, multiqc_config, [], multiqc_logo
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Count number of samples with the same samplename
def count_samples(ch_samples) {
    ch_samples.map { meta, fastq ->
        return [meta.samplename, [meta, fastq]]
    }
    // this should group per samplename
    .groupTuple()
    // Count the number of samples per samplename
    .map { samplename, meta_fastq ->
        count = meta_fastq.size()
        return [meta_fastq,count]
    }
    // split the meta_fastq list into channel items
    .transpose()
    // add the count variable to meta
    .map{ meta_fastq, count ->
        return [meta_fastq[0] + [count:count], meta_fastq[1]]
    }
}

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
    // dirty fix to have the `single_end` key in the meta map
    meta.single_end = false

    return [meta, bam ? bam : cram]
}

// Merge fastq meta with sample info
def merge_sample_info(ch_fastq, ch_sample_info) {
    ch_fastq
    .combine(ch_sample_info)
    .map { meta1, fastq, meta2 ->
        def meta = meta1.findAll{true}
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

        rg.ID = [fcid,lane].join(".")
        rg.PU = [fcid, lane, index].findAll().join(".")
        rg.PL = "ILLUMINA"
    } else if (fields.size() == 5) {
        fcid = fields[0]
        rg.ID = fcid
    }
    return rg
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
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
