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
            ? [file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true), file(params.multiqc_config, checkIfExists: true)]
            : [file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)],
        params.multiqc_logo ? file(params.multiqc_logo, checkIfExists: true) : [],
        params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
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
        params.hook_url,
        PREPROCESSING.out.multiqc_report,
    )

    publish:
    demultiplex_reports        = PREPROCESSING.out.demultiplex_reports.map { meta, report -> [meta, file(report).resolve()] }.transpose()
    demultiplex_logs           = PREPROCESSING.out.demultiplex_logs.map { meta, log -> [meta, file(log).resolve()] }.transpose()
    demultiplex_fastq          = PREPROCESSING.out.demultiplex_fastq.transpose()
    falco_html                 = PREPROCESSING.out.falco_html
    falco_txt                  = PREPROCESSING.out.falco_txt
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
    demultiplex_fastq {
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
