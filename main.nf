#!/usr/bin/env nextflow

// enable module feature
nextflow.enable.dsl = 2

// run workflow in workflows/
workflow {

    include { MULTI_EXPOSURE_MR } from './workflows/multi_exposure_mr'
    MULTI_EXPOSURE_MR()
}
