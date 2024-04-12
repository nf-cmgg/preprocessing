#!/usr/bin/env nextflow

// MODULES
include { MOSDEPTH              } from "../../../modules/nf-core/mosdepth/main.nf"
include { SAMTOOLS_COVERAGE     } from "../../../modules/nf-core/samtools/coverage/main"
include { PANELCOVERAGE         } from "../../../modules/local/panelcoverage/main"

workflow COVERAGE {
    take:
        ch_meta_cram_crai_fasta_fai_roi  // channel: [mandatory] [meta, cram, crai, fasta, fai, roi]
        ch_genelists                     // channel: [optional] [genelists]
    main:

        ch_versions         = Channel.empty()
        ch_coverageqc_files = Channel.empty()



        MOSDEPTH(
            ch_meta_cram_crai_fasta_fai_roi.map{
                meta, cram, crai, fasta, fai, roi ->
                    return [meta, cram, crai, roi, fasta]
                }
            )
        ch_versions = ch_versions.mix(MOSDEPTH.out.versions)

        SAMTOOLS_COVERAGE(ch_meta_cram_crai_fasta_fai_roi)
        ch_versions = ch_versions.mix(SAMTOOLS_COVERAGE.out.versions)
        ch_coverageqc_files = ch_coverageqc_files.merge(SAMTOOLS_COVERAGE.out.coverage)

        PANELCOVERAGE(
            MOSDEPTH.out.per_base_bed.join(MOSDEPTH.out.per_base_csi),
            ch_genelists
        )
        ch_versions = ch_versions.mix(PANELCOVERAGE.out.versions)
        ch_coverageqc_files = ch_coverageqc_files.mix(PANELCOVERAGE.out.regiondist)

    emit:
        mosdepth_summary  = MOSDEPTH.out.summary_txt
        mosdepth_global   = MOSDEPTH.out.global_txt
        mosdepth_regions  = MOSDEPTH.out.regions_txt
        samtools_coverage = SAMTOOLS_COVERAGE.out.coverage
        panelcoverage     = PANELCOVERAGE.out.regiondist
        versions          = ch_versions
}
