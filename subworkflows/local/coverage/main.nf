#!/usr/bin/env nextflow

// MODULES
include { MOSDEPTH          } from "../../../modules/nf-core/mosdepth/main.nf"
include { SAMTOOLS_COVERAGE } from "../../../modules/nf-core/samtools/coverage/main"
include { PANELCOVERAGE     } from "../../../modules/local/panelcoverage/main"

workflow COVERAGE {
    take:
    ch_meta_cram_crai_fasta_fai_roi // channel: [mandatory] [meta, cram, crai, fasta, fai, roi]
    ch_genelists // channel: [optional] [genelists]

    main:

    ch_versions = channel.empty()
    ch_coverageqc_files = channel.empty()

    MOSDEPTH(
        ch_meta_cram_crai_fasta_fai_roi.map { meta, cram, crai, fasta, _fai, roi ->
            return [meta, cram, crai, roi, fasta]
        }
    )

    SAMTOOLS_COVERAGE(
        ch_meta_cram_crai_fasta_fai_roi.map { meta, cram, crai, fasta, fai, _roi ->
            return [meta, cram, crai, fasta, fai]
        }
    )
    ch_coverageqc_files = ch_coverageqc_files.merge(SAMTOOLS_COVERAGE.out.coverage)

    PANELCOVERAGE(
        MOSDEPTH.out.per_base_bed.join(MOSDEPTH.out.per_base_csi).combine(ch_genelists).map { meta, bed, index, genelists ->
            // Because groovy typing sucks ass; apparently an array of 1 is automatically converted to a string...
            def genelists_array = genelists !instanceof List ? [genelists] : genelists
            def filtered_genelists = meta.tag.toLowerCase() == "seqcap"
                ? genelists_array.findAll { genelist -> genelist.name.toLowerCase().contains("seqcap") }
                : genelists_array.findAll { genelist -> !genelist.name.toLowerCase().contains("seqcap") }

            if (filtered_genelists.size() > 0) {
                return [
                    meta,
                    bed,
                    index,
                    filtered_genelists,
                ]
            }
        }
    )
    ch_coverageqc_files = ch_coverageqc_files.mix(PANELCOVERAGE.out.regiondist)

    emit:
    mosdepth_global         = MOSDEPTH.out.global_txt
    mosdepth_summary        = MOSDEPTH.out.summary_txt
    mosdepth_regions        = MOSDEPTH.out.regions_txt
    mosdepth_per_base_d4    = MOSDEPTH.out.per_base_d4
    mosdepth_per_base_bed   = MOSDEPTH.out.per_base_bed
    mosdepth_per_base_csi   = MOSDEPTH.out.per_base_csi
    mosdepth_regions_bed    = MOSDEPTH.out.regions_bed
    mosdepth_regions_csi    = MOSDEPTH.out.regions_csi
    mosdepth_quantized_bed  = MOSDEPTH.out.quantized_bed
    mosdepth_quantized_csi  = MOSDEPTH.out.quantized_csi
    mosdepth_thresholds_bed = MOSDEPTH.out.thresholds_bed
    mosdepth_thresholds_csi = MOSDEPTH.out.thresholds_csi
    samtools_coverage       = SAMTOOLS_COVERAGE.out.coverage
    panelcoverage           = PANELCOVERAGE.out.regiondist
}
