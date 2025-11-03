#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-cmgg/preprocessing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-cmgg/preprocessing
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

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
    demultiplex_interop = PREPROCESSING.out.demultiplex_interop
    demultiplex_reports = PREPROCESSING.out.demultiplex_reports.map { meta, reports -> [ meta, files("${reports.toUri()}/*") ] }.transpose(by:1)
    demultiplex_logs    = PREPROCESSING.out.demultiplex_logs.map { meta, logs -> [ meta, files("${logs.toUri()}/*") ] }.transpose(by:1)
}

output {
    // TODO also add the RunInfo.xml file as output, needs a module update
    demultiplex_interop { path "InterOp" }
    demultiplex_reports { path { meta, report ->
        report >> (meta.lane ? "Reports/LOO${meta.lane}/${report.name}" as String : "Reports/${report.name}")
    } }
    demultiplex_logs { path { meta, log ->
        log >> (meta.lane ? "Logs/LOO${meta.lane}/${log.name}" as String : "Logs/${log.name}")
    } }
}
