/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Modules
include { FASTP                  } from '../modules/nf-core/fastp/main'
include { MD5SUM                 } from '../modules/nf-core/md5sum/main'
include { MOSDEPTH               } from '../modules/nf-core/mosdepth/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { SAMTOOLS_COVERAGE      } from '../modules/nf-core/samtools/coverage/main'

// Subworkflows
include { BAM_QC                 } from '../subworkflows/local/bam_qc/main'
include { BCL_DEMULTIPLEX        } from '../subworkflows/nf-core/bcl_demultiplex/main'
include { COVERAGE               } from '../subworkflows/local/coverage/main'
include { FASTQ_TO_UCRAM         } from '../subworkflows/local/fastq_to_unaligned_cram/main'
include { FASTQ_TO_CRAM          } from '../subworkflows/local/fastq_to_aligned_cram/main'

// Functions
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_preprocessing_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PREPROCESSING {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    genomes        // map: genome reference files
    aligner        // string: global aligner to use
    markdup        // string: markdup method to use
    roi            // file: regions of interest bed file to be applied to all samples
    genelists      // file: directory containing genelist bed files for coverage analysis

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ch_samplesheet
    .branch {meta, fastq_1, fastq_2, samplesheet, sampleinfo, flowcell ->
        illumina_flowcell   : (flowcell && samplesheet && sampleinfo) && !(fastq_1 || fastq_2)
            return [meta, samplesheet, sampleinfo, flowcell]
        fastq               : (fastq_1) && !(flowcell || samplesheet || sampleinfo)
            return [meta, [fastq_1, fastq_2].findAll()]
        other: true
            error "Unable to determine input type, please check inputs"
    }
    .set{ch_inputs_from_samplesheet}

    roi = roi ? file(roi, checkIfExists:true) : null

    genelists = genelists ? file(genelists + "/*.bed", checkIfExists:true) : []

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// PROCESS FLOWCELL INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    ch_inputs_from_samplesheet.illumina_flowcell
    .multiMap { meta, samplesheet, sampleinfo, flowcell ->
        flowcell: [meta, samplesheet, flowcell]
        info    : sampleinfo
    }
    .set{ ch_illumina_flowcell }

    // BCL_DEMULTIPLEX([meta, samplesheet, flowcell], demultiplexer)
    BCL_DEMULTIPLEX(ch_illumina_flowcell.flowcell, "bclconvert")
    BCL_DEMULTIPLEX.out.fastq.dump(tag: "DEMULTIPLEX: fastq",pretty: true)
    ch_multiqc_files = ch_multiqc_files.mix(
        BCL_DEMULTIPLEX.out.reports.map { meta, reports -> return reports},
        BCL_DEMULTIPLEX.out.stats.map   { meta, stats   -> return stats  }
    )
    ch_versions = ch_versions.mix(BCL_DEMULTIPLEX.out.versions)

    // Add metadata to demultiplexed fastq's
    merge_sample_info(
        BCL_DEMULTIPLEX.out.fastq,
        parse_sample_info_csv(ch_illumina_flowcell.info)
    )
    .set {ch_demultiplexed_fastq}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// PROCESS FASTQ INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    ch_inputs_from_samplesheet.fastq
    .map { meta, fastq ->
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
    .set {ch_input_fastq}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ASSOCIATE CORRECT GENOME AND COUNT SAMPLE REPLICATES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
    ch_input_fastq
    .mix(ch_demultiplexed_fastq)
    // set genome based on organism key
    .map{ meta, reads ->
        if (meta.organism && !meta.genome) {
            switch (meta.organism) {
                case ~/(?i)Homo[\s_]sapiens/:
                    meta = meta + ["genome":"GRCh38"]
                    break
                // case ~/(?i)Mus[\s_]musculus/:
                //     meta = meta + ["genome":"GRCm39"]
                //     break
                // case ~/(?i)Danio[\s_]rerio/:
                //     meta = meta + ["genome":"GRCz11"]
                //     break
                default:
                    meta = meta + ["genome": null ]
                    break
            }
        }
        if (genomes && genomes[meta.genome]){
            meta = meta + ["genome": genomes[meta.genome]]
        }
        // set the aligner
        if (aligner && !meta.aligner) {
            meta = meta + ["aligner": aligner]
        }
        // set the ROI
        if (roi && !meta.roi) {
            meta = meta + ["roi": roi]
        }
        return [meta, reads]
    }
    // Count the number of samples per samplename
    .map{ meta, reads -> [meta.samplename, [meta, reads]]}
    .groupTuple()
    .map{ samplename, meta_fastq -> [meta_fastq, meta_fastq.size()]}
    .transpose()
    .map{meta_fastq, count -> [meta_fastq[0] + ['count': count], meta_fastq[1]]}
    // Clean up metadata
    .map{meta, fastq ->
        clean_meta = meta - meta.subMap('fcid','lane','library')
        return [clean_meta, fastq]
    }
    .set{ch_fastq_per_sample}

    ch_fastq_per_sample.dump(tag:"FASTQ per sample", pretty: true)


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// FASTQ TRIMMING AND QC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    // MODULE: fastp
    // Run QC, trimming and adapter removal
    // FASTP([meta, fastq], adapter_fasta, save_trimmed, save_merged)
    FASTP(ch_fastq_per_sample, [], false, false)
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.map { meta, json -> return json} )
    ch_versions      = ch_versions.mix(FASTP.out.versions)

    // edit meta.id to match sample name
    FASTP.out.reads
    .map { meta, reads ->
        def read_files = meta.single_end.toBoolean() ? reads : reads.sort{ a,b -> a.getName().tokenize('.')[0] <=> b.getName().tokenize('.')[0] }.collate(2)
        return [
            meta + [ chunks: read_files instanceof List ? read_files.size() : [read_files].size() ],
            read_files
        ]
    }
    // transpose to get read pairs
    .transpose()
    // set new meta.id to include split number
    .map { meta, reads ->
        def new_id = reads instanceof List ? reads[0].getName() - ~/_1.fastp.*/ : reads.getName() - ~/.fastp.*/
        return [
            meta - meta.subMap('id') + [ id: new_id ],
            reads
        ]
    }
    // split samples into human and non human data
    .branch { meta, reads ->
        supported: meta.genome
        other: true
    }
    .set { ch_trimmed_reads }

    ch_trimmed_reads.supported.dump(tag:"Supported trimmed reads per sample", pretty: true)
    ch_trimmed_reads.other.dump(tag:"Other trimmed reads per sample", pretty: true)


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// STEP: FASTQ TO UNALIGNED CRAM CONVERSION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    FASTQ_TO_UCRAM(ch_trimmed_reads.other)
    ch_versions = ch_versions.mix(FASTQ_TO_UCRAM.out.versions)

/*

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// STEP: FASTQ TO ALIGNED CRAM CONVERSION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
    ch_trimmed_reads.supported.map{ meta, reads ->
        return [
            meta,
            reads,
            meta.aligner,
            GenomeUtils.getGenomeAttribute(meta.genome, meta.aligner),
            GenomeUtils.getGenomeAttribute(meta.genome, "fasta")
            ]
    }
    .set{ch_meta_reads_aligner_index_fasta}

    FASTQ_TO_CRAM(
        ch_meta_reads_aligner_index_fasta,
        markdup
    )

    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_TO_CRAM.out.multiqc_files)
    ch_versions = ch_versions.mix(FASTQ_TO_CRAM.out.versions)


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// STEP: COVERAGE ANALYSIS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
    FASTQ_TO_CRAM.out.cram_crai
    .map { meta, cram, crai ->
        if (meta.roi) {
            return [
                meta,
                cram,
                crai,
                GenomeUtils.getGenomeAttribute(meta.genome, "fasta"),
                GenomeUtils.getGenomeAttribute(meta.genome, "fai"),
                file(meta.roi, checkIfExists:true),
            ]
        } else {
            return [
                meta,
                cram,
                crai,
                GenomeUtils.getGenomeAttribute(meta.genome, "fasta"),
                GenomeUtils.getGenomeAttribute(meta.genome, "fai"),
                [],
            ]
        }
    }
    .set{ch_cram_crai_fasta_fai_roi}

    COVERAGE(ch_cram_crai_fasta_fai_roi, genelists)
    ch_multiqc_files = ch_multiqc_files.mix(
        COVERAGE.out.mosdepth_summary   .map{ meta, txt -> return txt },
        COVERAGE.out.mosdepth_global    .map{ meta, txt -> return txt },
        COVERAGE.out.mosdepth_regions   .map{ meta, txt -> return txt },
        COVERAGE.out.samtools_coverage  .map{ meta, txt -> return txt },
    )
    ch_versions = ch_versions.mix(COVERAGE.out.versions)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// STEP: QC FOR ALIGNMENTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
    FASTQ_TO_CRAM.out.cram_crai
    .map { meta, cram, crai ->
        if (meta.roi) {
            return [
                meta,
                cram,
                crai,
                file(meta.roi, checkIfExists:true),
                GenomeUtils.getGenomeAttribute(meta.genome, "fasta"),
                GenomeUtils.getGenomeAttribute(meta.genome, "fai"),
                GenomeUtils.getGenomeAttribute(meta.genome, "dict"),
            ]
        } else {
            return [
                meta,
                cram,
                crai,
                [],
                GenomeUtils.getGenomeAttribute(meta.genome, "fasta"),
                GenomeUtils.getGenomeAttribute(meta.genome, "fai"),
                GenomeUtils.getGenomeAttribute(meta.genome, "dict"),
            ]
        }
    }
    .set{ch_cram_crai_roi_fasta_fai_dict}

    BAM_QC(ch_cram_crai_roi_fasta_fai_dict, params.disable_picard_metrics)
    ch_multiqc_files = ch_multiqc_files.mix(
        BAM_QC.out.samtools_stats           .map{ meta, txt -> return txt },
        BAM_QC.out.samtools_flagstat        .map{ meta, txt -> return txt },
        BAM_QC.out.samtools_idxstats        .map{ meta, txt -> return txt },
        BAM_QC.out.picard_multiplemetrics   .map{ meta, txt -> return txt },
        BAM_QC.out.picard_wgsmetrics        .map{ meta, txt -> return txt },
        BAM_QC.out.picard_hsmetrics         .map{ meta, txt -> return txt },
    )
    ch_versions = ch_versions.mix(BAM_QC.out.versions)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// STEP: CHECKSUMS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    MD5SUM(FASTQ_TO_CRAM.out.cram_crai.map{ meta, cram, crai -> return [meta,cram] })

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// AGGREGATE QC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
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
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
