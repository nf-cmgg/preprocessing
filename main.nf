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
        params.aligner,
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
        PREPROCESSING.out.multiqc_report,
    )

    publish:
    demultiplex_interop = PREPROCESSING.out.demultiplex_interop.transpose(by:1)
    demultiplex_reports = PREPROCESSING.out.demultiplex_reports.transpose(by:1)
    demultiplex_logs    = PREPROCESSING.out.demultiplex_logs.transpose(by:1)
    fastp_json          = PREPROCESSING.out.fastp_json
    fastp_html          = PREPROCESSING.out.fastp_html
    ucrams              = PREPROCESSING.out.ucrams
    crams               = PREPROCESSING.out.crams
    align_reports       = PREPROCESSING.out.align_reports
    sormadup_metrics    = PREPROCESSING.out.sormadup_metrics
    mosdepth_global = PREPROCESSING.out.mosdepth_global
    mosdepth_summary = PREPROCESSING.out.mosdepth_summary
    mosdepth_regions = PREPROCESSING.out.mosdepth_regions
    mosdepth_per_base_d4 = PREPROCESSING.out.mosdepth_per_base_d4
    mosdepth_per_base_bed = PREPROCESSING.out.mosdepth_per_base_bed
    mosdepth_per_base_csi = PREPROCESSING.out.mosdepth_per_base_csi
    mosdepth_regions_bed = PREPROCESSING.out.mosdepth_regions_bed
    mosdepth_regions_csi = PREPROCESSING.out.mosdepth_regions_csi
    mosdepth_quantized_bed = PREPROCESSING.out.mosdepth_quantized_bed
    mosdepth_quantized_csi = PREPROCESSING.out.mosdepth_quantized_csi
    mosdepth_thresholds_bed = PREPROCESSING.out.mosdepth_thresholds_bed
    mosdepth_thresholds_csi = PREPROCESSING.out.mosdepth_thresholds_csi
    samtools_coverage = PREPROCESSING.out.samtools_coverage
    panelcoverage = PREPROCESSING.out.panelcoverage
    samtools_stats = PREPROCESSING.out.samtools_stats
    samtools_flagstat = PREPROCESSING.out.samtools_flagstat
    samtools_idxstats = PREPROCESSING.out.samtools_idxstats
    picard_multiplemetrics = PREPROCESSING.out.picard_multiplemetrics
    picard_multiplemetrics_pdf = PREPROCESSING.out.picard_multiplemetrics_pdf
    picard_wgsmetrics = PREPROCESSING.out.picard_wgsmetrics
    picard_hsmetrics = PREPROCESSING.out.picard_hsmetrics
    md5sums = PREPROCESSING.out.md5sums
    multiqc_report = PREPROCESSING.out.multiqc_report
    multiqc_data = PREPROCESSING.out.multiqc_data
    multiqc_plots = PREPROCESSING.out.multiqc_plots

}

output {
    demultiplex_interop { path { _meta, bin ->
        bin >> "Interop/${bin.name}"
    } }
    demultiplex_reports { path { meta, report ->
        def out_path = meta.lane ? "Reports/L00${meta.lane}/${report.name}" as String : "Reports/${report.name}"
        report >> out_path
    } }
    demultiplex_logs { path { meta, log ->
        def out_path = meta.lane ? "Logs/L00${meta.lane}/${log.name}" as String : "Logs/${log.name}"
        log >> out_path
    } }
    fastp_json { path { meta, json ->
        json >> "${meta.samplename}/${json.name}"
    } }
    fastp_html { path { meta, html ->
        html >> "${meta.samplename}/${html.name}"
    } }
    ucrams { path { meta, cram ->
        cram >> "${meta.samplename}/${meta.samplename}.unaligned.cram"
    } }
    crams { path { meta, cram, crai ->
        cram >> "${meta.samplename}/${meta.samplename}.cram"
        crai >> "${meta.samplename}/${meta.samplename}.cram.crai"
    } }
    align_reports { path { meta, log ->
        log >> "${meta.samplename}/${log.name}"
    } }
    sormadup_metrics { path { meta, metrics ->
        metrics >> "${meta.samplename}/${meta.samplename}.duplicate_metrics.txt"
    } }
    mosdepth_global { path { meta, _file -> "${meta.samplename}/" } }
    mosdepth_summary { path { meta, _file -> "${meta.samplename}/" } }
    mosdepth_regions { path { meta, _file -> "${meta.samplename}/" } }
    mosdepth_per_base_d4 { path { meta, _file -> "${meta.samplename}/" } }
    mosdepth_per_base_bed { path { meta, _file -> "${meta.samplename}/" } }
    mosdepth_per_base_csi { path { meta, _file -> "${meta.samplename}/" } }
    mosdepth_regions_bed { path { meta, _file -> "${meta.samplename}/" } }
    mosdepth_regions_csi { path { meta, _file -> "${meta.samplename}/" } }
    mosdepth_quantized_bed { path { meta, _file -> "${meta.samplename}/" } }
    mosdepth_quantized_csi { path { meta, _file -> "${meta.samplename}/" } }
    mosdepth_thresholds_bed { path { meta, _file -> "${meta.samplename}/" } }
    mosdepth_thresholds_csi { path { meta, _file -> "${meta.samplename}/" } }
    samtools_coverage { path { meta, _file -> "${meta.samplename}/" } }
    panelcoverage { path { meta, _file -> "${meta.samplename}/" } }
    samtools_stats { path { meta, _file -> "${meta.samplename}/" } }
    samtools_flagstat { path { meta, _file -> "${meta.samplename}/" } }
    samtools_idxstats { path { meta, _file -> "${meta.samplename}/" } }
    picard_multiplemetrics { path { meta, _file -> "${meta.samplename}/" } }
    picard_multiplemetrics_pdf { path { meta, _file -> "${meta.samplename}/" } }
    picard_wgsmetrics { path { meta, _file -> "${meta.samplename}/" } }
    picard_hsmetrics { path { meta, _file -> "${meta.samplename}/" } }
    md5sums { path { meta, _file -> "${meta.samplename}/" } }
    multiqc_report { path "multiqc/" }
    multiqc_data { path "multiqc/" }
    multiqc_plots { path "multiqc/" }
}
