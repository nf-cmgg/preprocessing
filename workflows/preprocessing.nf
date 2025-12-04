include { samplesheetToList      } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Modules
include { FASTP                         } from '../modules/nf-core/fastp/main'
include { MD5SUM                        } from '../modules/nf-core/md5sum/main'
include { MOSDEPTH                      } from '../modules/nf-core/mosdepth/main'
include { MULTIQC as MULTIQC_LIBRARY    } from '../modules/nf-core/multiqc/main'
include { MULTIQC as MULTIQC_MAIN       } from '../modules/nf-core/multiqc/main'
include { SAMTOOLS_COVERAGE             } from '../modules/nf-core/samtools/coverage/main'

// Subworkflows
include { BAM_QC                        } from '../subworkflows/local/bam_qc/main'
include { BCL_DEMULTIPLEX               } from '../subworkflows/nf-core/bcl_demultiplex/main'
include { COVERAGE                      } from '../subworkflows/local/coverage/main'
include { FASTQ_TO_UCRAM                } from '../subworkflows/local/fastq_to_unaligned_cram/main'
include { FASTQ_TO_CRAM                 } from '../subworkflows/local/fastq_to_aligned_cram/main'

// Functions
include { paramsSummaryMap              } from 'plugin/nf-schema'
include { paramsSummaryMultiqc          } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText        } from '../subworkflows/local/utils_nfcore_preprocessing_pipeline'
include { getGenomeAttribute            } from '../subworkflows/local/utils_nfcore_preprocessing_pipeline'

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
    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    ch_samplesheet
        .branch { meta, fastq_1, fastq_2, samplesheet, sampleinfo, flowcell ->
            illumina_flowcell: (flowcell && samplesheet && sampleinfo) && !(fastq_1 || fastq_2)
            return [meta, samplesheet, sampleinfo, flowcell]
            fastq: (fastq_1) && !(flowcell || samplesheet || sampleinfo)
            return [meta, [fastq_1, fastq_2].findAll()]
            other: true
            error("Unable to determine input type, please check inputs")
        }
        .set { ch_inputs_from_samplesheet }

    roi = roi ? file(roi, checkIfExists: true) : null

    genelists = genelists ? channel.value(file(genelists + "/*.bed", checkIfExists: true)) : channel.empty()

    /*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// PROCESS FLOWCELL INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    ch_inputs_from_samplesheet.illumina_flowcell
        .multiMap { meta, samplesheet, sampleinfo, flowcell ->
            flowcell: [meta, samplesheet, flowcell]
            info: samplesheetToList(sampleinfo, "assets/schema_sampleinfo.json")
        }
        .set { ch_illumina_flowcell }

    // BCL_DEMULTIPLEX([meta, samplesheet, flowcell], demultiplexer)
    BCL_DEMULTIPLEX(ch_illumina_flowcell.flowcell, "bclconvert")
    BCL_DEMULTIPLEX.out.fastq.dump(tag: "DEMULTIPLEX: fastq", pretty: true)
    ch_multiqc_files = ch_multiqc_files.mix(
        BCL_DEMULTIPLEX.out.reports,
        BCL_DEMULTIPLEX.out.stats
    )
    ch_versions = ch_versions.mix(BCL_DEMULTIPLEX.out.versions)

    BCL_DEMULTIPLEX.out.fastq
        .map { meta, fastq -> [meta.samplename, meta, fastq] }
        .set { ch_demultiplexed_fastq }

    ch_illumina_flowcell.info
        .flatten()
        .transpose()
        .map { sampleinfo -> [sampleinfo.samplename, sampleinfo] }
        .set { ch_sampleinfo }

    // Merge fastq meta with sample info
    ch_demultiplexed_fastq
        .combine(ch_sampleinfo, by: 0)
        .map { samplename, meta, fastq, sampleinfo ->
            def new_meta = meta + sampleinfo
            def readgroup = readgroup_from_fastq(fastq[0])
            readgroup = readgroup + ['SM': samplename, 'LB': new_meta.library ?: ""]
            new_meta = new_meta + ['readgroup': readgroup]
            return [new_meta, fastq]
        }
        .groupTuple(by: [0])
        .map { meta, fq ->
            return [meta, fq.flatten().unique()]
        }
        .set { ch_demultiplexed_fastq_with_sampleinfo }
    /*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// PROCESS FASTQ INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    ch_inputs_from_samplesheet.fastq
        .map { meta, fastq ->
            // if no fastq_2, then single-end
            def single_end = fastq[1] ? false : true
            // add readgroup metadata
            def rg = readgroup_from_fastq(fastq[0])
            // if the sample name starts with "snp_", remove it so the sampletracking works later on.
            def samplename = meta.samplename.startsWith("snp_") ? meta.samplename.substring(4) : meta.samplename
            rg = rg + [
                'SM': samplename,
                'LB': meta.library ?: "",
                'PL': meta.platform ?: rg.PL,
                'ID': meta.readgroup ?: rg.ID,
            ]
            def meta_with_readgroup = meta + ['single_end': single_end, 'readgroup': rg]
            return [meta_with_readgroup, fastq]
        }
        .set { ch_input_fastq }

    /*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ASSOCIATE CORRECT GENOME AND COUNT SAMPLE REPLICATES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
    ch_input_fastq
        .mix(ch_demultiplexed_fastq_with_sampleinfo)
        .map { meta, reads ->
            if (meta.organism && !meta.genome) {
                if (meta.organism ==~ /(?i)Homo[\s_]sapiens/) {
                    meta = meta + ["genome": "GRCh38"]
                }
                else if (meta.organism ==~ /(?i)Mus[\s_]musculus/) {
                    meta = meta + ["genome": "mm10"]
                }
                else if (meta.organism ==~ /(?i)Danio[\s_]rerio/) {
                    meta = meta + ["genome": "GRCz11"]
                }
                else {
                    meta = meta + ["genome": null]
                }
            }
            if (genomes && genomes[meta.genome]) {
                meta = meta + ["genome_data": genomes[meta.genome]]
            }
            else {
                meta = meta + ["genome_data": [:]]
            }
            // set the aligner
            if (aligner && !meta.aligner) {
                meta = meta + ["aligner": aligner]
            }
            // set the ROI
            // // Special case for coPGT samples
            // // if there's no global ROI AND no sample speficic ROI
            // // AND the sample tag is "coPGT-M", set the sample ROI to "roi_copgt"
            if (!roi && !meta.roi && meta.tag == "coPGT-M") {
                meta = meta + ["roi": getGenomeAttribute(meta.genome_data, "roi_copgt")]
            }
            // // if there's a global ROI AND no sample specific ROI
            // // set the global ROI to the sample
            if (roi && !meta.roi) {
                meta = meta + ["roi": roi]
            }
            return [meta, reads]
        }
        .map { meta, reads -> [meta.samplename, [meta, reads]] }
        .groupTuple()
        .map { _samplename, meta_fastq -> [meta_fastq, meta_fastq.size()] }
        .transpose()
        .map { meta_fastq, count -> [meta_fastq[0] + ['count': count], meta_fastq[1]] }
        .map { meta, fastq ->
            return [meta - meta.subMap('fcid', 'lane'), fastq]
        }
        .set { ch_fastq_per_sample }

    ch_fastq_per_sample.dump(tag: "FASTQ per sample", pretty: true)


    /*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// FASTQ TRIMMING AND QC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    // MODULE: fastp
    // Run QC, trimming and adapter removal
    // FASTP([meta, fastq], adapter_fasta, save_trimmed, save_merged)
    FASTP(ch_fastq_per_sample.map{ meta, fastq -> return [meta, fastq, []] }, false, false, false)
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json)
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    // edit meta.id to match sample name
    FASTP.out.reads
        .map { meta, reads ->
            def read_files = meta.single_end.toBoolean() ? reads : reads.sort { a, b -> a.getName().tokenize('.')[0] <=> b.getName().tokenize('.')[0] }.collate(2)
            return [
                meta + [chunks: read_files instanceof List ? read_files.size() : [read_files].size()],
                read_files,
            ]
        }
        .transpose()
        .map { meta, reads ->
            def new_id = reads instanceof List ? reads[0].getName() - ~/_1.fastp.*/ : reads.getName() - ~/.fastp.*/
            return [
                meta - meta.subMap('id') + [id: new_id],
                reads,
            ]
        }
        .branch { meta, _reads ->
            supported: meta.genome_data instanceof Map && meta.genome_data.size() > 0
            other: true
        }
        .set { ch_trimmed_reads }

    ch_trimmed_reads.supported.dump(tag: "Supported trimmed reads per sample", pretty: true)
    ch_trimmed_reads.other.dump(tag: "Other trimmed reads per sample", pretty: true)


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
    ch_trimmed_reads.supported
        .map { meta, reads ->
            return [
                meta,
                reads,
                meta.aligner,
                getGenomeAttribute(meta.genome_data, meta.aligner),
                getGenomeAttribute(meta.genome_data, "fasta"),
                getGenomeAttribute(meta.genome_data, "gtf"),
            ]
        }
        .set { ch_meta_reads_aligner_index_fasta_gtf }

    FASTQ_TO_CRAM(
        ch_meta_reads_aligner_index_fasta_gtf,
        markdup,
    )

    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_TO_CRAM.out.sormadup_metrics)
    ch_versions = ch_versions.mix(FASTQ_TO_CRAM.out.versions)


    /*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// STEP: FILTER SAMPLES WITH 'SNP' TAG
// samples with SNP tag contain only data for sample tracking
// and as such don't need all the QC steps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    FASTQ_TO_CRAM.out.cram_crai
        .filter { meta, _cram, _crai ->
            meta.tag != "SNP"
        }
        .set { ch_no_snp_samples }

    /*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// STEP: COVERAGE ANALYSIS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
    ch_no_snp_samples
        .map { meta, cram, crai ->
            if (meta.roi) {
                return [
                    meta,
                    cram,
                    crai,
                    getGenomeAttribute(meta.genome_data, "fasta"),
                    getGenomeAttribute(meta.genome_data, "fai"),
                    file(meta.roi, checkIfExists: true),
                ]
            }
            else {
                return [
                    meta,
                    cram,
                    crai,
                    getGenomeAttribute(meta.genome_data, "fasta"),
                    getGenomeAttribute(meta.genome_data, "fai"),
                    [],
                ]
            }
        }
        .set { ch_cram_crai_fasta_fai_roi }

    def mosdepth_global_out = channel.empty()
    def mosdepth_summary_out = channel.empty()
    def mosdepth_regions_out = channel.empty()
    def mosdepth_per_base_d4_out = channel.empty()
    def mosdepth_per_base_bed_out = channel.empty()
    def mosdepth_per_base_csi_out = channel.empty()
    def mosdepth_regions_bed_out = channel.empty()
    def mosdepth_regions_csi_out = channel.empty()
    def mosdepth_quantized_bed_out = channel.empty()
    def mosdepth_quantized_csi_out = channel.empty()
    def mosdepth_thresholds_bed_out = channel.empty()
    def mosdepth_thresholds_csi_out = channel.empty()
    def samtools_coverage_out = channel.empty()
    def panelcoverage_out = channel.empty()
    if (params.run_coverage) {
        COVERAGE(ch_cram_crai_fasta_fai_roi, genelists)
        ch_multiqc_files = ch_multiqc_files.mix(
            COVERAGE.out.mosdepth_summary,
            COVERAGE.out.mosdepth_global,
            COVERAGE.out.mosdepth_regions,
            COVERAGE.out.samtools_coverage,
        )
        mosdepth_global_out = COVERAGE.out.mosdepth_global
        mosdepth_summary_out = COVERAGE.out.mosdepth_summary
        mosdepth_regions_out = COVERAGE.out.mosdepth_regions
        mosdepth_per_base_d4_out = COVERAGE.out.mosdepth_per_base_d4
        mosdepth_per_base_bed_out = COVERAGE.out.mosdepth_per_base_bed
        mosdepth_per_base_csi_out = COVERAGE.out.mosdepth_per_base_csi
        mosdepth_regions_bed_out = COVERAGE.out.mosdepth_regions_bed
        mosdepth_regions_csi_out = COVERAGE.out.mosdepth_regions_csi
        mosdepth_quantized_bed_out = COVERAGE.out.mosdepth_quantized_bed
        mosdepth_quantized_csi_out = COVERAGE.out.mosdepth_quantized_csi
        mosdepth_thresholds_bed_out = COVERAGE.out.mosdepth_thresholds_bed
        mosdepth_thresholds_csi_out = COVERAGE.out.mosdepth_thresholds_csi
        samtools_coverage_out = COVERAGE.out.samtools_coverage
        panelcoverage_out = COVERAGE.out.panelcoverage
        ch_versions = ch_versions.mix(COVERAGE.out.versions)
    }

    /*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// STEP: QC FOR ALIGNMENTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
    ch_no_snp_samples
        .map { meta, cram, crai ->
            if (meta.roi) {
                return [
                    meta,
                    cram,
                    crai,
                    file(meta.roi, checkIfExists: true),
                    getGenomeAttribute(meta.genome_data, "fasta"),
                    getGenomeAttribute(meta.genome_data, "fai"),
                    getGenomeAttribute(meta.genome_data, "dict"),
                ]
            }
            else {
                return [
                    meta,
                    cram,
                    crai,
                    [],
                    getGenomeAttribute(meta.genome_data, "fasta"),
                    getGenomeAttribute(meta.genome_data, "fai"),
                    getGenomeAttribute(meta.genome_data, "dict"),
                ]
            }
        }
        .set { ch_cram_crai_roi_fasta_fai_dict }

    BAM_QC(ch_cram_crai_roi_fasta_fai_dict, params.disable_picard_metrics)
    ch_multiqc_files = ch_multiqc_files.mix(
        BAM_QC.out.samtools_stats,
        BAM_QC.out.samtools_flagstat,
        BAM_QC.out.samtools_idxstats,
        BAM_QC.out.picard_multiplemetrics,
        BAM_QC.out.picard_wgsmetrics,
        BAM_QC.out.picard_wgsmetrics,
        BAM_QC.out.picard_hsmetrics,
    )
    ch_versions = ch_versions.mix(BAM_QC.out.versions)

    /*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// STEP: CHECKSUMS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    MD5SUM(
        FASTQ_TO_CRAM.out.cram_crai.map { meta, cram, _crai ->
            return [meta, cram]
        },
        false,
    )
    ch_versions = ch_versions.mix(MD5SUM.out.versions.first())

    /*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// AGGREGATE QC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
    //
    // Collate and save software versions
    //
    def topic_versions = Channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_cmgg_preprocessing_software_mqc_versions.yml',
            sort: true,
            newLine: true
        )
        .map { file -> [[id: 'main'], file] } // add meta for multiqc
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ? channel.fromPath(params.multiqc_config, checkIfExists: true) : channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ? channel.fromPath(params.multiqc_logo, checkIfExists: true) : channel.empty()

    summary_params           = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml').map{ file -> [[id: 'main'], file] })

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true).map{ file -> [[id: 'main'], file] })

    ch_multiqc_files = ch_multiqc_files.map { meta, files ->
        return [meta.library ? [id: meta.library] : [id: 'main'], files]
    }
    .groupTuple()
    .branch { meta, files ->
        main: meta.id == 'main'
            return files.flatten()
        library: meta.id != 'main'
            return [meta, files.flatten()]
    }
    ch_multiqc_files.main.dump(tag: "MULTIQC files - main", pretty: true)
    ch_multiqc_files.library.dump(tag: "MULTIQC files - library", pretty: true)

    MULTIQC_MAIN(
        ch_multiqc_files.main.collect().map{ files -> [[id: 'main'], files] },
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        [],
    )

    MULTIQC_LIBRARY(
        ch_multiqc_files.library,
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        [],
    )

    emit:
    demultiplex_interop         = BCL_DEMULTIPLEX.out.interop
    demultiplex_reports         = BCL_DEMULTIPLEX.out.reports
    demultiplex_logs            = BCL_DEMULTIPLEX.out.logs
    fastp_json                  = FASTP.out.json
    fastp_html                  = FASTP.out.html
    ucrams                      = FASTQ_TO_UCRAM.out.cram
    crams                       = FASTQ_TO_CRAM.out.cram_crai
    align_reports               = FASTQ_TO_CRAM.out.align_reports
    sormadup_metrics            = FASTQ_TO_CRAM.out.sormadup_metrics
    mosdepth_global             = mosdepth_global_out
    mosdepth_summary            = mosdepth_summary_out
    mosdepth_regions            = mosdepth_regions_out
    mosdepth_per_base_d4        = mosdepth_per_base_d4_out
    mosdepth_per_base_bed       = mosdepth_per_base_bed_out
    mosdepth_per_base_csi       = mosdepth_per_base_csi_out
    mosdepth_regions_bed        = mosdepth_regions_bed_out
    mosdepth_regions_csi        = mosdepth_regions_csi_out
    mosdepth_quantized_bed      = mosdepth_quantized_bed_out
    mosdepth_quantized_csi      = mosdepth_quantized_csi_out
    mosdepth_thresholds_bed     = mosdepth_thresholds_bed_out
    mosdepth_thresholds_csi     = mosdepth_thresholds_csi_out
    samtools_coverage           = samtools_coverage_out
    panelcoverage               = panelcoverage_out
    samtools_stats              = BAM_QC.out.samtools_stats
    samtools_flagstat           = BAM_QC.out.samtools_flagstat
    samtools_idxstats           = BAM_QC.out.samtools_idxstats
    picard_multiplemetrics      = BAM_QC.out.picard_multiplemetrics
    picard_multiplemetrics_pdf  = BAM_QC.out.picard_multiplemetrics_pdf
    picard_wgsmetrics           = BAM_QC.out.picard_wgsmetrics
    picard_hsmetrics            = BAM_QC.out.picard_hsmetrics
    md5sums                     = MD5SUM.out.checksum
    multiqc_main_report         = MULTIQC_MAIN.out.report.toList()
    multiqc_main_data           = MULTIQC_MAIN.out.data
    multiqc_main_plots          = MULTIQC_MAIN.out.plots
    multiqc_library_report      = MULTIQC_LIBRARY.out.report.toList()
    multiqc_library_data        = MULTIQC_LIBRARY.out.data
    multiqc_library_plots       = MULTIQC_LIBRARY.out.plots
    versions                    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// https://github.com/nf-core/sarek/blob/7ba61bde8e4f3b1932118993c766ed33b5da465e/workflows/sarek.nf#L1014-L1040
def readgroup_from_fastq(path) {
    // expected format:
    // xx:yy:FLOWCELLID:LANE:... (seven fields)
    // or
    // FLOWCELLID:LANE:xx:... (five fields)
    def line

    path.withInputStream { fq ->
        def gzipStream = new java.util.zip.GZIPInputStream(fq) as InputStream
        def decoder = new InputStreamReader(gzipStream, 'ASCII')
        def buffered = new BufferedReader(decoder)
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
        // def sequencer_serial = fields[0]
        // def run_number       = fields[1]
        def fcid = fields[2]
        def lane = fields[3]
        def index = fields[-1] =~ /[GATC+-]/ ? fields[-1] : ""

        rg.ID = [fcid, lane].join(".")
        rg.PU = [fcid, lane, index].findAll().join(".")
        rg.PL = "ILLUMINA"
    }
    else if (fields.size() == 5) {
        def fcid = fields[0]
        rg.ID = fcid
    }
    return rg
}
