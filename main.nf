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
        params.markdup,
        params.roi,
        params.genelists,
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
        PREPROCESSING.out.multiqc_main_report,
    )

    publish:
    demultiplex_interop        = PREPROCESSING.out.demultiplex_interop.transpose(by: 1)
    demultiplex_reports        = PREPROCESSING.out.demultiplex_reports.transpose(by: 1)
    demultiplex_logs           = PREPROCESSING.out.demultiplex_logs.transpose(by: 1)
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
    multiqc_main_report        = PREPROCESSING.out.multiqc_main_report
    multiqc_main_data          = PREPROCESSING.out.multiqc_main_data
    multiqc_main_plots         = PREPROCESSING.out.multiqc_main_plots
    multiqc_library_report     = PREPROCESSING.out.multiqc_library_report
    multiqc_library_data       = PREPROCESSING.out.multiqc_library_data
    multiqc_library_plots      = PREPROCESSING.out.multiqc_library_plots
}

output {
    demultiplex_interop {
        path { _meta, bin ->
            bin >> "Interop/${bin.name}"
        }
    }
    demultiplex_reports {
        path { meta, report ->
            def out_path = meta.lane ? "Reports/L00${meta.lane}/${report.name}" as String : "Reports/${report.name}"
            report >> out_path
        }
    }
    demultiplex_logs {
        path { meta, log ->
            def out_path = meta.lane ? "Logs/L00${meta.lane}/${log.name}" as String : "Logs/${log.name}"
            log >> out_path
        }
    }
    demultiplex_fastq {
        path { meta, fastq ->
            def out_path = meta.library ? "${meta.library}/${meta.samplename}/${fastq.name}" as String : "${meta.samplename}/${fastq.name}"
            fastq >> out_path
        }
    }
    falco_html {
        path { meta, html ->
            def out_path = meta.library ? "${meta.library}/${meta.samplename}/${html.name}" as String : "${meta.samplename}/${html.name}"
            html >> out_path
        }
    }
    falco_txt {
        path { meta, txt ->
            def out_path = meta.library ? "${meta.library}/${meta.samplename}/${txt.name}" as String : "${meta.samplename}/${txt.name}"
            txt >> out_path
        }
    }
    fastp_json {
        path { meta, json ->
            def out_path = meta.library ? "${meta.library}/${meta.samplename}/${json.name}" as String : "${meta.samplename}/${json.name}"
            json >> out_path
        }
    }
    fastp_html {
        path { meta, html ->
            def out_path = meta.library ? "${meta.library}/${meta.samplename}/${html.name}" as String : "${meta.samplename}/${html.name}"
            html >> out_path
        }
    }
    crams {
        path { meta, cram, crai ->
            def out_cram = meta.library ? "${meta.library}/${meta.samplename}/${meta.samplename}.cram" as String : "${meta.samplename}/${meta.samplename}.cram"
            def out_crai = meta.library ? "${meta.library}/${meta.samplename}/${meta.samplename}.cram.crai" as String : "${meta.samplename}/${meta.samplename}.cram.crai"
            cram >> out_cram
            crai >> out_crai
        }
    }
    rna_splice_junctions {
        path { meta, sjt ->
            def out_path = meta.library ? "${meta.library}/${meta.samplename}/${sjt.name}" as String : "${meta.samplename}/${sjt.name}"
            sjt >> out_path
        }
    }
    rna_junctions {
        path { meta, junctions ->
            def out_path = meta.library ? "${meta.library}/${meta.samplename}/${junctions.name}" as String : "${meta.samplename}/${junctions.name}"
            junctions >> out_path
        }
    }
    align_reports {
        path { meta, log ->
            def out_path = meta.library ? "${meta.library}/${meta.samplename}/${log.name}" as String : "${meta.samplename}/${log.name}"
            log >> out_path
        }
    }
    sormadup_metrics {
        path { meta, metrics ->
            def out_path = meta.library ? "${meta.library}/${meta.samplename}/${meta.samplename}.duplicate_metrics.txt" as String : "${meta.samplename}/${meta.samplename}.duplicate_metrics.txt"
            metrics >> out_path
        }
    }
    mosdepth_global {
        path { meta, _file ->
            def out_path = meta.library ? "${meta.library}/${meta.samplename}/" as String : "${meta.samplename}/"
            return out_path
        }
    }
    mosdepth_summary {
        path { meta, _file ->
            def out_path = meta.library ? "${meta.library}/${meta.samplename}/" as String : "${meta.samplename}/"
            return out_path
        }
    }
    mosdepth_regions {
        path { meta, _file ->
            def out_path = meta.library ? "${meta.library}/${meta.samplename}/" as String : "${meta.samplename}/"
            return out_path
        }
    }
    mosdepth_per_base_d4 {
        path { meta, _file ->
            def out_path = meta.library ? "${meta.library}/${meta.samplename}/" as String : "${meta.samplename}/"
            return out_path
        }
    }
    mosdepth_per_base_bed {
        path { meta, _file ->
            def out_path = meta.library ? "${meta.library}/${meta.samplename}/" as String : "${meta.samplename}/"
            return out_path
        }
    }
    mosdepth_per_base_csi {
        path { meta, _file ->
            def out_path = meta.library ? "${meta.library}/${meta.samplename}/" as String : "${meta.samplename}/"
            return out_path
        }
    }
    mosdepth_regions_bed {
        path { meta, _file ->
            def out_path = meta.library ? "${meta.library}/${meta.samplename}/" as String : "${meta.samplename}/"
            return out_path
        }
    }
    mosdepth_regions_csi {
        path { meta, _file ->
            def out_path = meta.library ? "${meta.library}/${meta.samplename}/" as String : "${meta.samplename}/"
            return out_path
        }
    }
    mosdepth_quantized_bed {
        path { meta, _file ->
            def out_path = meta.library ? "${meta.library}/${meta.samplename}/" as String : "${meta.samplename}/"
            return out_path
        }
    }
    mosdepth_quantized_csi {
        path { meta, _file ->
            def out_path = meta.library ? "${meta.library}/${meta.samplename}/" as String : "${meta.samplename}/"
            return out_path
        }
    }
    mosdepth_thresholds_bed {
        path { meta, _file ->
            def out_path = meta.library ? "${meta.library}/${meta.samplename}/" as String : "${meta.samplename}/"
            return out_path
        }
    }
    mosdepth_thresholds_csi {
        path { meta, _file ->
            def out_path = meta.library ? "${meta.library}/${meta.samplename}/" as String : "${meta.samplename}/"
            return out_path
        }
    }
    samtools_coverage {
        path { meta, _file ->
            def out_path = meta.library ? "${meta.library}/${meta.samplename}/" as String : "${meta.samplename}/"
            return out_path
        }
    }
    panelcoverage {
        path { meta, _file ->
            def out_path = meta.library ? "${meta.library}/${meta.samplename}/" as String : "${meta.samplename}/"
            return out_path
        }
    }
    samtools_stats {
        path { meta, _file ->
            def out_path = meta.library ? "${meta.library}/${meta.samplename}/" as String : "${meta.samplename}/"
            return out_path
        }
    }
    samtools_flagstat {
        path { meta, _file ->
            def out_path = meta.library ? "${meta.library}/${meta.samplename}/" as String : "${meta.samplename}/"
            return out_path
        }
    }
    samtools_idxstats {
        path { meta, _file ->
            def out_path = meta.library ? "${meta.library}/${meta.samplename}/" as String : "${meta.samplename}/"
            return out_path
        }
    }
    picard_multiplemetrics {
        path { meta, _file ->
            def out_path = meta.library ? "${meta.library}/${meta.samplename}/" as String : "${meta.samplename}/"
            return out_path
        }
    }
    picard_multiplemetrics_pdf {
        path { meta, _file ->
            def out_path = meta.library ? "${meta.library}/${meta.samplename}/" as String : "${meta.samplename}/"
            return out_path
        }
    }
    picard_wgsmetrics {
        path { meta, _file ->
            def out_path = meta.library ? "${meta.library}/${meta.samplename}/" as String : "${meta.samplename}/"
            return out_path
        }
    }
    picard_hsmetrics {
        path { meta, _file ->
            def out_path = meta.library ? "${meta.library}/${meta.samplename}/" as String : "${meta.samplename}/"
            return out_path
        }
    }
    md5sums {
        path { meta, _file ->
            def out_path = meta.library ? "${meta.library}/${meta.samplename}/" as String : "${meta.samplename}/"
            return out_path
        }
    }
    multiqc_main_report {
        path "multiqc/"
    }
    multiqc_main_data {
        path "multiqc/"
    }
    multiqc_main_plots {
        path "multiqc/"
    }
    multiqc_library_report {
        path "multiqc/"
    }
    multiqc_library_data {
        path "multiqc/"
    }
    multiqc_library_plots {
        path "multiqc/"
    }
}
