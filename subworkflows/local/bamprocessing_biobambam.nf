#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.options = [:]

include {BIOBAMBAM_BAMCAT} from "../../modules/nf-core/modules/biobambam/bamcat"
include {BIOBAMBAM_BAMSORMADUP} from "../../modules/nf-core/modules/biobambam/bamsormadup"
include {BIOBAMBAM_BAMMERGE} from "../../modules/nf-core/modules/biobambam/bammerge"
