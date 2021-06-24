#!/usr/bin/env nextflow

// enable module feature
nextflow.enable.dsl = 2

include { MULTI_EXPOSURE_MR } from './workflows/multi_exposure_mr'

// run workflow in workflows/
workflow {
    MULTI_EXPOSURE_MR()
}
