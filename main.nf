#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-cmgg/preprocessing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-cmgg/preprocessing
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_preprocessing_pipeline'
include { PREPROCESSING           } from './workflows/preprocessing'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_preprocessing_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params {

    // Path to comma-separated or yaml file containing information about the samples in the experiment.
    input: Path

    // The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.
    outdir: Path

    // Email address for completion summary.
    email: String?

    // MultiQC report title. Printed as page header, used for filename if not otherwise specified.
    multiqc_title: String?

    genomes: Map = [:]

    // Specify how many reads each split of a FastQ file contains. Set 0 to turn off splitting at all.
    split_fastq: Integer = 100000000

    // Directory containing gene list bed files for granular coverage analysis
    genelists: Path?

    // Git commit id for Institutional configs.
    custom_config_version: String = 'main'

    // Base directory for custom configs.
    custom_config_base: String = 'https://raw.githubusercontent.com/nf-cmgg/configs/main'

    // Institutional config name.
    config_profile_name: String?

    // Institutional config description.
    config_profile_description: String?

    // Institutional config contact information.
    config_profile_contact: String?

    // Institutional config URL link.
    config_profile_url: String?

    // Display version and exit.
    version: Boolean = false

    // Method used to save pipeline results to output directory.
    publish_dir_mode: String = 'copy'

    // Email address for completion summary, only when pipeline fails.
    email_on_fail: String?

    // Send plain-text email instead of HTML.
    plaintext_email: Boolean = false

    // File size limit when attaching MultiQC reports to summary emails.
    max_multiqc_email_size: String = '25.MB'

    // Do not use coloured log outputs.
    monochrome_logs: Boolean = false

    // Incoming hook URL for messaging service
    hook_url: String = System.getenv('HOOK_URL')

    // Custom config file to supply to MultiQC.
    multiqc_config: Path?

    // Custom logo file to supply to MultiQC. File name must also be set in the MultiQC config file
    multiqc_logo: Path?

    // Custom MultiQC yaml file containing HTML including a methods description.
    multiqc_methods_description: Path?

    // Boolean whether to validate parameters against the schema at runtime
    validate_params: Boolean = true

    // Base URL or local path to location of pipeline test dataset files
    pipelines_testdata_base_path: String = 'https://raw.githubusercontent.com/nf-core/test-datasets/'

    // Suffix to add to the trace report filename. Default is the date and time in the format yyyy-MM-dd_HH-mm-ss.
    trace_report_suffix: String

    // Display the help message.
    help = false

    // Display the full detailed help message.
    help_full: Boolean = false

    // Display hidden parameters in the help message (only works when --help or --help_full are provided).
    show_hidden: Boolean = false

    // Directory / URL base for iGenomes references.
    igenomes_base: String = '/references/'

    // Do not load the iGenomes reference config.
    igenomes_ignore: Boolean = false

    // Name of iGenomes reference.
    genome: String?
}

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION(
        params.version,
        params.validate_params,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden,
    )

    //
    // WORKFLOW: Run main workflow
    //
    PREPROCESSING(
        PIPELINE_INITIALISATION.out.samplesheet,
        params.genomes,
        params.genelists,
        params.multiqc_config
            ? [file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true), params.multiqc_config]
            : [file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)],
        params.multiqc_logo ? params.multiqc_logo : [],
        params.multiqc_methods_description ? params.multiqc_methods_description : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true),
        params.outdir,
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        PREPROCESSING.out.multiqc_report,
    )

    publish:
    demultiplex_reports        = PREPROCESSING.out.demultiplex_reports.transpose()
    demultiplex_logs           = PREPROCESSING.out.demultiplex_logs.transpose()
    demultiplex_interop        = PREPROCESSING.out.demultiplex_interop.transpose(by: 1)
    fastq                      = PREPROCESSING.out.fastq.transpose()
    falco_html                 = PREPROCESSING.out.falco_html.transpose()
    falco_txt                  = PREPROCESSING.out.falco_txt.transpose()
    fastp_json                 = PREPROCESSING.out.fastp_json
    fastp_html                 = PREPROCESSING.out.fastp_html
    crams                      = PREPROCESSING.out.crams
    rna_splice_junctions       = PREPROCESSING.out.rna_splice_junctions
    rna_junctions              = PREPROCESSING.out.rna_junctions
    align_reports              = PREPROCESSING.out.align_reports
    sormadup_metrics           = PREPROCESSING.out.sormadup_metrics
    mosdepth_global            = PREPROCESSING.out.mosdepth_global
    mosdepth_summary           = PREPROCESSING.out.mosdepth_summary
    mosdepth_regions           = PREPROCESSING.out.mosdepth_regions
    mosdepth_per_base_d4       = PREPROCESSING.out.mosdepth_per_base_d4
    mosdepth_per_base_bed      = PREPROCESSING.out.mosdepth_per_base_bed
    mosdepth_per_base_csi      = PREPROCESSING.out.mosdepth_per_base_csi
    mosdepth_regions_bed       = PREPROCESSING.out.mosdepth_regions_bed
    mosdepth_regions_csi       = PREPROCESSING.out.mosdepth_regions_csi
    mosdepth_quantized_bed     = PREPROCESSING.out.mosdepth_quantized_bed
    mosdepth_quantized_csi     = PREPROCESSING.out.mosdepth_quantized_csi
    mosdepth_thresholds_bed    = PREPROCESSING.out.mosdepth_thresholds_bed
    mosdepth_thresholds_csi    = PREPROCESSING.out.mosdepth_thresholds_csi
    samtools_coverage          = PREPROCESSING.out.samtools_coverage
    panelcoverage              = PREPROCESSING.out.panelcoverage
    samtools_stats             = PREPROCESSING.out.samtools_stats
    samtools_flagstat          = PREPROCESSING.out.samtools_flagstat
    samtools_idxstats          = PREPROCESSING.out.samtools_idxstats
    picard_multiplemetrics     = PREPROCESSING.out.picard_multiplemetrics
    picard_multiplemetrics_pdf = PREPROCESSING.out.picard_multiplemetrics_pdf
    picard_wgsmetrics          = PREPROCESSING.out.picard_wgsmetrics
    picard_hsmetrics           = PREPROCESSING.out.picard_hsmetrics
    md5sums                    = PREPROCESSING.out.md5sums
    multiqc_report             = PREPROCESSING.out.multiqc_report
    multiqc_data               = PREPROCESSING.out.multiqc_data
    multiqc_plots              = PREPROCESSING.out.multiqc_plots
    multiqcsav_report          = PREPROCESSING.out.multiqcsav_report
    multiqcsav_data            = PREPROCESSING.out.multiqcsav_data
    multiqcsav_plots           = PREPROCESSING.out.multiqcsav_plots
}

output {
    demultiplex_reports {
        path { meta, report ->
            report >> (meta.lane ? "Reports/L00${meta.lane}/${report.name}" : "Reports/${report.name}")
        }
    }
    demultiplex_logs {
        path { meta, log ->
            log >> (meta.lane ? "Logs/L00${meta.lane}/${log.name}" : "Logs/${log.name}")
        }
    }
    demultiplex_interop {
        path { _meta, bin ->
            bin >> "Interop/${bin.name}"
        }
    }
    fastq {
        path { meta, fastq ->
            fastq >> (meta.library ? "${meta.library}/${meta.samplename}/${fastq.name}" : "${meta.samplename}/${fastq.name}")
        }
    }
    falco_html {
        path { meta, html ->
            html >> (meta.library ? "${meta.library}/${meta.samplename}/${html.name}" : "${meta.samplename}/${html.name}")
        }
    }
    falco_txt {
        path { meta, txt ->
            txt >> (meta.library ? "${meta.library}/${meta.samplename}/${txt.name}" : "${meta.samplename}/${txt.name}")
        }
    }
    fastp_json {
        path { meta, json ->
            json >> (meta.library ? "${meta.library}/${meta.samplename}/${json.name}" : "${meta.samplename}/${json.name}")
        }
    }
    fastp_html {
        path { meta, html ->
            html >> (meta.library ? "${meta.library}/${meta.samplename}/${html.name}" : "${meta.samplename}/${html.name}")
        }
    }
    crams {
        path { meta, cram, crai ->
            cram >> (meta.library ? "${meta.library}/${meta.samplename}/${meta.samplename}.cram" : "${meta.samplename}/${meta.samplename}.cram")
            crai >> (meta.library ? "${meta.library}/${meta.samplename}/${meta.samplename}.cram.crai" : "${meta.samplename}/${meta.samplename}.cram.crai")
        }
    }
    rna_splice_junctions {
        path { meta, sjt ->
            sjt >> (meta.library ? "${meta.library}/${meta.samplename}/${sjt.name}" : "${meta.samplename}/${sjt.name}")
        }
    }
    rna_junctions {
        path { meta, junctions ->
            junctions >> (meta.library ? "${meta.library}/${meta.samplename}/${junctions.name}" : "${meta.samplename}/${junctions.name}")
        }
    }
    align_reports {
        path { meta, log ->
            log >> (meta.library ? "${meta.library}/${meta.samplename}/${log.name}" : "${meta.samplename}/${log.name}")
        }
    }
    sormadup_metrics {
        path { meta, metrics ->
            metrics >> (meta.library ? "${meta.library}/${meta.samplename}/${meta.samplename}.duplicate_metrics.txt" : "${meta.samplename}/${meta.samplename}.duplicate_metrics.txt")
        }
    }
    mosdepth_global {
        path { meta, _file ->
            return (meta.library ? "${meta.library}/${meta.samplename}/" : "${meta.samplename}/")
        }
    }
    mosdepth_summary {
        path { meta, _file ->
            return (meta.library ? "${meta.library}/${meta.samplename}/" : "${meta.samplename}/")
        }
    }
    mosdepth_regions {
        path { meta, _file ->
            return (meta.library ? "${meta.library}/${meta.samplename}/" : "${meta.samplename}/")
        }
    }
    mosdepth_per_base_d4 {
        path { meta, _file ->
            return (meta.library ? "${meta.library}/${meta.samplename}/" : "${meta.samplename}/")
        }
    }
    mosdepth_per_base_bed {
        path { meta, _file ->
            return (meta.library ? "${meta.library}/${meta.samplename}/" : "${meta.samplename}/")
        }
    }
    mosdepth_per_base_csi {
        path { meta, _file ->
            return (meta.library ? "${meta.library}/${meta.samplename}/" : "${meta.samplename}/")
        }
    }
    mosdepth_regions_bed {
        path { meta, _file ->
            return (meta.library ? "${meta.library}/${meta.samplename}/" : "${meta.samplename}/")
        }
    }
    mosdepth_regions_csi {
        path { meta, _file ->
            return (meta.library ? "${meta.library}/${meta.samplename}/" : "${meta.samplename}/")
        }
    }
    mosdepth_quantized_bed {
        path { meta, _file ->
            return (meta.library ? "${meta.library}/${meta.samplename}/" : "${meta.samplename}/")
        }
    }
    mosdepth_quantized_csi {
        path { meta, _file ->
            return (meta.library ? "${meta.library}/${meta.samplename}/" : "${meta.samplename}/")
        }
    }
    mosdepth_thresholds_bed {
        path { meta, _file ->
            return (meta.library ? "${meta.library}/${meta.samplename}/" : "${meta.samplename}/")
        }
    }
    mosdepth_thresholds_csi {
        path { meta, _file ->
            return (meta.library ? "${meta.library}/${meta.samplename}/" : "${meta.samplename}/")
        }
    }
    samtools_coverage {
        path { meta, _file ->
            return (meta.library ? "${meta.library}/${meta.samplename}/" : "${meta.samplename}/")
        }
    }
    panelcoverage {
        path { meta, _file ->
            return (meta.library ? "${meta.library}/${meta.samplename}/" : "${meta.samplename}/")
        }
    }
    samtools_stats {
        path { meta, _file ->
            return (meta.library ? "${meta.library}/${meta.samplename}/" : "${meta.samplename}/")
        }
    }
    samtools_flagstat {
        path { meta, _file ->
            return (meta.library ? "${meta.library}/${meta.samplename}/" : "${meta.samplename}/")
        }
    }
    samtools_idxstats {
        path { meta, _file ->
            return (meta.library ? "${meta.library}/${meta.samplename}/" : "${meta.samplename}/")
        }
    }
    picard_multiplemetrics {
        path { meta, _file ->
            return (meta.library ? "${meta.library}/${meta.samplename}/" : "${meta.samplename}/")
        }
    }
    picard_multiplemetrics_pdf {
        path { meta, _file ->
            return (meta.library ? "${meta.library}/${meta.samplename}/" : "${meta.samplename}/")
        }
    }
    picard_wgsmetrics {
        path { meta, _file ->
            return (meta.library ? "${meta.library}/${meta.samplename}/" : "${meta.samplename}/")
        }
    }
    picard_hsmetrics {
        path { meta, _file ->
            return (meta.library ? "${meta.library}/${meta.samplename}/" : "${meta.samplename}/")
        }
    }
    md5sums {
        path { meta, _file ->
            return (meta.library ? "${meta.library}/${meta.samplename}/" : "${meta.samplename}/")
        }
    }
    multiqcsav_report {
        path "multiqc/"
    }
    multiqcsav_data {
        path "multiqc/"
    }
    multiqcsav_plots {
        path "multiqc/"
    }
    multiqc_report {
        path { meta, _file ->
            return (meta.id ? "${meta.id}/multiqc/" : "multiqc/")
        }
    }
    multiqc_data {
        path { meta, _file ->
            return (meta.id ? "${meta.id}/multiqc/" : "multiqc/")
        }
    }
    multiqc_plots {
        path { meta, _file ->
            return (meta.id ? "${meta.id}/multiqc/" : "multiqc/")
        }
    }
}
