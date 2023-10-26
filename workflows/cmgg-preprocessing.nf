#!/usr/bin/env nextflow

import nextflow.Nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

// Validate input parameters
WorkflowCmggpreprocessing.initialise(params, log)

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
include { BCL_DEMULTIPLEX   } from "../subworkflows/nf-core/bcl_demultiplex/main"
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
include { MOSDEPTH                    } from "../modules/nf-core/mosdepth/main"
include { MULTIQC                     } from "../modules/nf-core/multiqc/main"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow CMGGPREPROCESSING {

    // input channels
    ch_input = params.input ? Channel.fromSamplesheet('input') : Channel.empty()

    // input values
    aligner = params.aligner
    markdup = params.markdup
    genomes = params.genomes

    // input options
    run_coverage   = params.run_coverage
    disable_picard = params.disable_picard_metrics

    // output channels
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // PROCESS INPUT CHANNEL
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    ch_input
    | branch {meta, bam, cram, fastq_1, fastq_2, samplesheet, sampleinfo, flowcell ->
        fc      : (flowcell && samplesheet && sampleinfo) && !(bam || cram || fastq_1 || fastq_2)
            return [meta, samplesheet, sampleinfo, flowcell]
        fastq   : (fastq_1) && !(bam || cram || flowcell || samplesheet || sampleinfo)
            return [meta, [fastq_1, fastq_2].findAll()]
        aligned : (bam || cram) && !(fastq_1 || fastq_2 || flowcell || samplesheet || sampleinfo)
            return [meta, [bam, cram].find()]
        other: true
            error "Unable to determine input type, please check inputs"
    }
    | set { ch_input }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// PROCESS FLOWCELL INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    ch_input.fc
    | multiMap { meta, samplesheet, sampleinfo, flowcell ->
        fc   : [meta, samplesheet, flowcell]
        info : sampleinfo
    }
    | set { ch_flowcell }

    // BCL_DEMULTIPLEX([meta, samplesheet, flowcell], demultiplexer)
    BCL_DEMULTIPLEX(ch_flowcell.fc, "bclconvert")
    BCL_DEMULTIPLEX.out.fastq.dump(tag: "DEMULTIPLEX: fastq",pretty: true)
    ch_multiqc_files = ch_multiqc_files.mix(
        BCL_DEMULTIPLEX.out.reports.map { meta, reports -> return reports},
        BCL_DEMULTIPLEX.out.stats.map   { meta, stats   -> return stats  }
    )
    ch_versions = ch_versions.mix(BCL_DEMULTIPLEX.out.versions)

    // Add metadata to demultiplexed fastq's
    merge_sample_info(
        BCL_DEMULTIPLEX.out.fastq,
        parse_sample_info_csv(ch_flowcell.info)
    )
    | set {ch_demultiplexed_fastq}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// PROCESS FASTQ INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    ch_input.fastq
    | map { meta, fastq ->
        // if no fastq_2, then single-end
        single_end = fastq[1] ? false : true
        // add readgroup metadata
        rg = readgroup_from_fastq(fastq[0])
        rg = rg + [ 'SM': meta.samplename,
                    'LB': meta.library ?: "",
                    'PL': meta.platform ?: rg.PL,
                    'ID': meta.readgroup ?: rg.ID
                ]

        meta_with_readgroup = meta + ['single_end': single_end, 'readgroup': rg]

        return [meta_with_readgroup, fastq]
    }
    | set {ch_input_fastq}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ASSOCIATE CORRECT GENOME TO SAMPLES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
    ch_input_fastq
    | mix(ch_demultiplexed_fastq, ch_input.aligned)
    | map { meta, reads ->
        // set genome based on organism key
        if (meta.organism && !meta.genome) {
            switch (meta.organism) {
                case ~/(?i)Homo[\s_]sapiens/:
                    meta = meta + ["genome":"GRCh38"]
                    break
                case ~/(?i)Mus[\s_]musculus/:
                    meta = meta + ["genome":"GRCm39"]
                    break
                case ~/(?i)Danio[\s_]rerio/:
                    meta = meta + ["genome":"GRCz11"]
                    break
                default:
                    meta = meta + ["genome": null ]
                    break
            }
        }
        if (genomes && genomes[meta.genome]){
            meta = meta + ["genome": genomes[meta.genome]]
        }
        return [meta, reads]
    }
    | branch {meta, reads ->
        fastq   : (reads instanceof List && reads[0] ==~ /.*.f(ast)?q\.gz$/) || (reads ==~ /.*.f(ast)?q\.gz$/)
        aligned : (reads instanceof List && reads[0] ==~ /.*.(b|cr)am$)/ )   || (reads ==~ /.*.(b|cr)am$/)
    }
    | set {ch_inputs}

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // STEP: PROCESS BAM/CRAM INPUTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Convert bam/cram inputs to fastq
    ch_inputs.aligned
    | map { meta, bam ->
        return [
            meta,
            bam,
            WorkflowMain.getGenomeAttribute(meta.genome, "fasta") ? WorkflowMain.getGenomeAttribute(meta.genome, "fasta") : [],
            WorkflowMain.getGenomeAttribute(meta.genome, "fai") ? WorkflowMain.getGenomeAttribute(meta.genome, "fai") : []
        ]
    }
    | set { ch_meta_reads_fasta_fai }

    BAM_TO_FASTQ(ch_meta_reads_fasta_fai)

    ch_converted_fastq = BAM_TO_FASTQ.out.fastq
    ch_converted_fastq.dump(tag: "MAIN: converted_fastq",pretty: true)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // STEP:GATHER PROCESSED INPUTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // "Gather" fastq's from demultiplex, cramtofastq and fastq inputs
    ch_inputs.fastq.mix(ch_converted_fastq)
    // count the number of samples per samplename
    | map { meta, fastq ->
        return [meta.samplename, [meta, fastq]]
    }
    // this should group per samplename
    | groupTuple()
    // Count the number of samples per samplename
    | map { samplename, meta_fastq ->
        count = meta_fastq.size()
        return [meta_fastq,count]
    }
    // split the meta_fastq list into channel items
    | transpose()
    // add the count variable to meta
    | map { meta_fastq, count ->
        return [meta_fastq[0] + [count:count], meta_fastq[1]]
    }
    | set { ch_sample_fastqs }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// FASTQ TRIMMING AND QC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    // MODULE: fastp
    // Run QC, trimming and adapter removal
    // FASTP([meta, fastq], adapter_fasta, save_trimmed, save_merged)
    FASTP(ch_sample_fastqs, [], false, false)
    FASTP.out.reads.dump(tag: "MAIN: fastp trimmed reads",pretty: true)
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.map { meta, json -> return json} )
    ch_versions      = ch_versions.mix(FASTP.out.versions)

    // edit meta.id to match sample name
    FASTP.out.reads
    | map { meta, reads ->
        def read_files = meta.single_end.toBoolean() ? reads : reads.sort{ a,b -> a.getName().tokenize('.')[0] <=> b.getName().tokenize('.')[0] }.collate(2)
        return [
            meta + [ chunks: read_files instanceof List ? read_files.size() : [read_files].size() ],
            read_files
        ]
    }
    // transpose to get read pairs
    | transpose()
    // set new meta.id to include split number
    | map { meta, reads ->
        def new_id = reads instanceof List ? reads[0].getName() - ~/_1.fastp.*/ : reads.getName() - ~/.fastp.*/
        return [
            meta - meta.subMap('id') + [ id: new_id ],
            reads
        ]
    }
    // split samples into human and non human data
    | branch { meta, reads ->
        supported: meta.genome
        other: true
    }
    | set { ch_trimmed_reads }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// STEP: FASTQ TO UNALIGNED CRAM CONVERSION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    FASTQ_TO_UCRAM(ch_trimmed_reads.other)
    ch_versions = ch_versions.mix(FASTQ_TO_UCRAM.out.versions)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// STEP: FASTQ TO ALIGNED CRAM CONVERSION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    FASTQ_TO_CRAM(
        ch_trimmed_reads.supported,
        aligner,
        markdup
    )

    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_TO_CRAM.out.multiqc_files)
    ch_versions = ch_versions.mix(FASTQ_TO_CRAM.out.versions)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// STEP: COVERAGE ANALYSIS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    // Generate coverage metrics and beds for each sample
    FASTQ_TO_CRAM.out.cram_crai
    | map { meta, cram, crai ->
        if (meta.roi) {
            return [
                meta,
                cram,
                crai,
                file(meta.roi, checkIfExists:true)
            ]
        } else {
            return [
                meta,
                cram,
                crai,
                []
            ]
        }
    }
    | set {ch_cram_crai_roi }

    if (run_coverage){
        ch_cram_crai_roi
        | multiMap { meta, cram, crai, roi ->
            meta_cram_crai_roi : [meta, cram, crai, roi]
            meta_fasta: [meta, WorkflowMain.getGenomeAttribute(meta.genome, "fasta")]
        }
        | set {ch_to_mosdepth}

        MOSDEPTH(
            ch_to_mosdepth.meta_cram_crai_roi,
            ch_to_mosdepth.meta_fasta
        )
        ch_multiqc_files = ch_multiqc_files.mix(
            MOSDEPTH.out.summary_txt,
            MOSDEPTH.out.global_txt,
            MOSDEPTH.out.regions_txt
        )
        ch_versions = ch_versions.mix(MOSDEPTH.out.versions).first()
    }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// STEP: QC FOR ALIGNMENTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    // Gather metrics from bam files
    // BAM_QC([meta, bam, bai], disable_picard)
    BAM_QC(FASTQ_TO_CRAM.out.cram_crai, disable_picard)
    ch_multiqc_files = ch_multiqc_files.mix( BAM_QC.out.metrics.map { meta, metrics -> return metrics} )
    ch_versions      = ch_versions.mix(BAM_QC.out.versions)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// STEP: REPORTING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    // MODULE: CUSTOM_DUMPSOFTWAREVERSIONS
    // Gather software versions for QC report
    CUSTOM_DUMPSOFTWAREVERSIONS (ch_versions.unique().collectFile(name: "collated_versions.yml"))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowCmggpreprocessing.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowCmggpreprocessing.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Parse sample info input map
def parse_sample_info_csv(csv_file) {
    csv_file.splitCsv(header: true, strip: true).map { row ->
        // check mandatory fields
        if (!(row.samplename)) log.error "Missing samplename field in sample info file"
        return row
    }
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
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
