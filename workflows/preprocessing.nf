include { samplesheetToList           } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Modules
include { BCLCONVERT                  } from '../modules/nf-core/bclconvert'
include { FALCO                       } from '../modules/nf-core/falco'
include { FASTP                       } from '../modules/nf-core/fastp'
include { MD5SUM                      } from '../modules/nf-core/md5sum'
include { MULTIQC                     } from '../modules/nf-core/multiqc'
include { MULTIQCSAV                  } from '../modules/nf-core/multiqcsav'
include { SAMTOOLS_COVERAGE           } from '../modules/nf-core/samtools/coverage'

// Subworkflows
include { BAM_QC                      } from '../subworkflows/local/bam_qc'
include { COVERAGE                    } from '../subworkflows/local/coverage'
include { FASTQ_TO_CRAM               } from '../subworkflows/local/fastq_to_aligned_cram'

// Functions
include { getReadgroupsFromBclconvert } from '../subworkflows/local/utils_nfcmgg_preprocessing_pipeline'
include { getReadgroupFromFastq       } from '../subworkflows/local/utils_nfcmgg_preprocessing_pipeline'
include { paramsSummaryMap            } from 'plugin/nf-schema'
include { paramsSummaryMultiqc        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText      } from '../subworkflows/local/utils_nfcore_preprocessing_pipeline'
include { getGenomeAttribute          } from '../subworkflows/local/utils_nfcore_preprocessing_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PREPROCESSING {
    take:
    ch_samplesheet // channel: samplesheet read in from --input
    genomes // map: genome reference files
    genelists // file: directory containing genelist bed files for coverage analysis
    multiqc_config // file(s): MultiQC config file(s)
    multiqc_logo // file: MultiQC logo file
    multiqc_methods_description // file: custom methods description for MultiQC report
    outdir // directory: output directory for the workflow results

    main:
    ch_multiqc_files = channel.empty()

    ch_samplesheet
        .branch { meta, fastq_1, fastq_2, samplesheet, sampleinfo, flowcell ->
            illumina_flowcell: (flowcell && samplesheet && sampleinfo) && !(fastq_1 || fastq_2)
            return [["id": meta.id, "lane": meta.lane], samplesheet, sampleinfo, flowcell]
            fastq: (fastq_1) && !(flowcell || samplesheet || sampleinfo)
            return [meta, [fastq_1, fastq_2].findAll()]
            other: true
            error("Unable to determine input type, please check inputs")
        }
        .set { ch_inputs_from_samplesheet }
    // construct a value channel containing an array of files, because the coverage subworkflow expects a channel of arrays of genelist files (to allow for multiple genelist files per sample)
    ch_genelists = genelists ? channel.fromPath(genelists + "/*.bed").collect().map { files -> [files] } : channel.empty()

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

    // BCLCONVERT([meta, samplesheet, flowcell])
    BCLCONVERT(ch_illumina_flowcell.flowcell)
    BCLCONVERT.out.fastq.dump(tag: "DEMULTIPLEX: fastq", pretty: true)

    getReadgroupsFromBclconvert(
        BCLCONVERT.out.reports.map { meta, reports ->
            return [meta, file(reports).resolve("fastq_list.csv")]
        },
        BCLCONVERT.out.fastq,
    ).dump(tag: "DEMULTIPLEX: fastq with meta", pretty: true).map { meta, fastq -> [meta.readgroup.SM, meta, fastq] }.set { ch_demultiplexed_fastq }

    // Run QC
    ch_mqcsav_input = ch_illumina_flowcell.flowcell
        .map { meta, _samplesheet, flowcell ->
            def interop = files(flowcell.resolve("InterOp/*.bin"), checkIfExists: true)
            def xml = files(flowcell.resolve("*.xml"), checkIfExists: true)
            return [meta, xml, interop]
        }
        .join(BCLCONVERT.out.reports, by: 0)
        .map { meta, xml, interop, reports ->
            return [meta - meta.subMap(['lane']), xml, interop, reports]
        }
        .groupTuple(by: [0])
        .map { meta, xml, interop, reports ->
            return [meta, xml.flatten().unique(), interop.flatten().unique(), reports.flatten(), multiqc_config, multiqc_logo, [], []]
        }
        .dump(tag: "MULTIQC SAV input", pretty: true)

    MULTIQCSAV(
        ch_mqcsav_input
    )

    // Merge fastq meta with sample info
    ch_illumina_flowcell.info
        .flatten()
        .transpose()
        .map { sampleinfo -> [sampleinfo.samplename, sampleinfo] }
        .set { ch_sampleinfo }

    ch_demultiplexed_fastq
        .combine(ch_sampleinfo, by: 0)
        .map { _samplename, meta, fastq, sampleinfo ->
            def new_rg = [:]
            if (sampleinfo.library) {
                new_rg = meta.readgroup + ['LB': sampleinfo.library]
            }
            else {
                new_rg = meta.readgroup
            }
            def new_meta = meta + sampleinfo + ['readgroup': new_rg]
            return [new_meta, fastq]
        }
        .groupTuple(by: [0])
        .map { meta, fq ->
            return [meta, fq.flatten().unique()]
        }
        .branch { meta, _fastq ->
            to_align: meta.aligner && meta.aligner != "false"
            other: true
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
            // if the sample name starts with "snp_", remove it so the sampletracking works later on.
            def samplename = meta.samplename.startsWith("snp_") ? meta.samplename.substring(4) : meta.samplename
            def rg = getReadgroupFromFastq(fastq[0], samplename, meta.library, meta.platform)
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
        .mix(ch_demultiplexed_fastq_with_sampleinfo.to_align)
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
        .branch { meta, _reads ->
            supported: meta.genome_data instanceof Map && meta.genome_data.size() > 0 && meta.aligner
            other: true
        }
        .set { ch_fastq_per_sample }

    ch_fastq_per_sample.supported.dump(tag: "Supported FASTQ per sample", pretty: true)
    ch_fastq_per_sample.other.dump(tag: "Other FASTQ per sample", pretty: true)

    /*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// FASTQ TRIMMING AND QC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    // MODULE: FALCO
    // Run FALCO for "unsupported" fastq QC
    // FALCO([meta, fastq])
    FALCO(ch_fastq_per_sample.other)
    ch_multiqc_files = ch_multiqc_files.mix(FALCO.out.html)
    ch_multiqc_files = ch_multiqc_files.mix(FALCO.out.txt)

    // MODULE: fastp
    // Run QC, trimming and adapter removal
    // FASTP([meta, fastq, adapter_fasta], save_trimmed, save_merged)
    FASTP(
        ch_fastq_per_sample.supported.map { meta, fastq ->
            return [meta, fastq, []]
        },
        false,
        false,
        false,
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json)

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
            def new_id = reads instanceof List ? reads[0].getName() - ~/_R1.fastp.*/ : reads.getName() - ~/.fastp.*/
            return [
                meta - meta.subMap('id') + [id: new_id],
                reads,
            ]
        }
        .set { ch_trimmed_reads }

    ch_trimmed_reads.dump(tag: "Supported trimmed reads per sample", pretty: true)

    /*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// STEP: FASTQ TO ALIGNED CRAM CONVERSION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
    ch_trimmed_reads
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
        ch_meta_reads_aligner_index_fasta_gtf
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_TO_CRAM.out.sormadup_metrics)

    /*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// STEP: COVERAGE ANALYSIS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
    FASTQ_TO_CRAM.out.cram_crai
        .filter { meta, _cram, _crai ->
            meta.run_coverage.toBoolean()
        }
        .map { meta, cram, crai ->
            return [
                meta,
                cram,
                crai,
                getGenomeAttribute(meta.genome_data, "fasta"),
                getGenomeAttribute(meta.genome_data, "fai"),
                meta.roi && meta.roi != [] ? file(meta.roi, checkIfExists: true) : [],
            ]
        }
        .set { ch_coverage }

    COVERAGE(ch_coverage, ch_genelists)
    ch_multiqc_files = ch_multiqc_files.mix(
        COVERAGE.out.mosdepth_summary,
        COVERAGE.out.mosdepth_global,
        COVERAGE.out.mosdepth_regions,
        COVERAGE.out.samtools_coverage,
    )

    /*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// STEP: QC FOR ALIGNMENTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
    FASTQ_TO_CRAM.out.cram_crai
        .map { meta, cram, crai ->
            return [
                meta,
                cram,
                crai,
                meta.roi && meta.roi != [] ? file(meta.roi, checkIfExists: true) : [],
                getGenomeAttribute(meta.genome_data, "fasta"),
                getGenomeAttribute(meta.genome_data, "fai"),
                getGenomeAttribute(meta.genome_data, "dict"),
            ]
        }
        .set { ch_bam_qc }

    BAM_QC(ch_bam_qc)
    ch_multiqc_files = ch_multiqc_files.mix(
        BAM_QC.out.samtools_stats,
        BAM_QC.out.samtools_flagstat,
        BAM_QC.out.samtools_idxstats,
        BAM_QC.out.picard_multiplemetrics,
        BAM_QC.out.picard_wgsmetrics,
        BAM_QC.out.picard_wgsmetrics,
        BAM_QC.out.picard_hsmetrics,
    )

    /*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// STEP: CHECKSUMS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    MD5SUM(
        ch_fastq_per_sample.other.mix(
            FASTQ_TO_CRAM.out.cram_crai.map { meta, cram, _crai ->
                return [meta, cram]
            }
        ),
        false,
    )

    /*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// AGGREGATE QC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [process[process.lastIndexOf(':') + 1..-1], "  ${tool}: ${version}"]
        }
        .groupTuple(by: 0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(topic_versions.versions_file)
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir.toUriString()}/pipeline_info",
            name: 'nf_cmgg_preprocessing_software_mqc_versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    // summary files without meta, e.g. versions, params
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_methods_description = channel.value(multiqc_methods_description ? methodsDescriptionText(multiqc_methods_description) : "")

    ch_summary_files = channel.empty()
        .mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        .mix(ch_collated_versions)
        .mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
        .toList()
        .map { files -> [files] }
        .dump(tag: "Summary files for MultiQC", pretty: true)

    ch_multiqc_input = ch_multiqc_files
        .map { meta, files ->
            def new_meta = meta.library ? [id: meta.library] : [id: 'multiqc']
            return [new_meta, files]
        }
        .groupTuple(by: 0)
        .combine(ch_summary_files)
        .map { meta, multiqc_files, summary_files ->
            return [meta, (multiqc_files + summary_files).flatten(), multiqc_config.flatten(), multiqc_logo, [], []]
        }
        .dump(tag: "MULTIQC files", pretty: true)


    // MULTIQC([meta, multiqc_files, multiqc_config, multiqc_logo, replace_names, sample_names])
    MULTIQC(ch_multiqc_input)

    emit:
    demultiplex_reports        = BCLCONVERT.out.reports.map { meta, reports ->
        return [meta, files(reports.resolve("*"))]
    }
    demultiplex_logs           = BCLCONVERT.out.logs.map { meta, logs ->
        return [meta, files(logs.resolve("*"))]
    }
    demultiplex_interop        = BCLCONVERT.out.interop
    demultiplex_fastq          = ch_demultiplexed_fastq_with_sampleinfo.other
    falco_html                 = FALCO.out.html
    falco_txt                  = FALCO.out.txt
    fastp_json                 = FASTP.out.json
    fastp_html                 = FASTP.out.html
    crams                      = FASTQ_TO_CRAM.out.cram_crai
    rna_splice_junctions       = FASTQ_TO_CRAM.out.rna_splice_junctions
    rna_junctions              = FASTQ_TO_CRAM.out.rna_junctions
    align_reports              = FASTQ_TO_CRAM.out.align_reports
    sormadup_metrics           = FASTQ_TO_CRAM.out.sormadup_metrics
    mosdepth_global            = COVERAGE.out.mosdepth_global
    mosdepth_summary           = COVERAGE.out.mosdepth_summary
    mosdepth_regions           = COVERAGE.out.mosdepth_regions
    mosdepth_per_base_d4       = COVERAGE.out.mosdepth_per_base_d4
    mosdepth_per_base_bed      = COVERAGE.out.mosdepth_per_base_bed
    mosdepth_per_base_csi      = COVERAGE.out.mosdepth_per_base_csi
    mosdepth_regions_bed       = COVERAGE.out.mosdepth_regions_bed
    mosdepth_regions_csi       = COVERAGE.out.mosdepth_regions_csi
    mosdepth_quantized_bed     = COVERAGE.out.mosdepth_quantized_bed
    mosdepth_quantized_csi     = COVERAGE.out.mosdepth_quantized_csi
    mosdepth_thresholds_bed    = COVERAGE.out.mosdepth_thresholds_bed
    mosdepth_thresholds_csi    = COVERAGE.out.mosdepth_thresholds_csi
    samtools_coverage          = COVERAGE.out.samtools_coverage
    panelcoverage              = COVERAGE.out.panelcoverage
    samtools_stats             = BAM_QC.out.samtools_stats
    samtools_flagstat          = BAM_QC.out.samtools_flagstat
    samtools_idxstats          = BAM_QC.out.samtools_idxstats
    picard_multiplemetrics     = BAM_QC.out.picard_multiplemetrics
    picard_multiplemetrics_pdf = BAM_QC.out.picard_multiplemetrics_pdf
    picard_wgsmetrics          = BAM_QC.out.picard_wgsmetrics
    picard_hsmetrics           = BAM_QC.out.picard_hsmetrics
    md5sums                    = MD5SUM.out.checksum
    multiqcsav_report          = MULTIQCSAV.out.report.toList()
    multiqcsav_data            = MULTIQCSAV.out.data.toList()
    multiqcsav_plots           = MULTIQCSAV.out.plots.toList()
    multiqc_report             = MULTIQC.out.report
    multiqc_data               = MULTIQC.out.data
    multiqc_plots              = MULTIQC.out.plots
}
