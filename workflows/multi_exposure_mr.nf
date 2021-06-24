
// import local modules
include { GET_INPUT_VARIANTS } from "$baseDir/modules/get_input_variants.nf"
include { CIS } from "$baseDir/workflows/workflow_split.nf"
include { ALL } from "$baseDir/workflows/workflow_split.nf"

// main pipeline
workflow MULTI_EXPOSURE_MR {

    // get exposure file names
    exposure_input_data = Channel.fromPath("${params.exposure.input_dir}/*")

    // reformat data into TwoSampleMR format w/ significant variants & clump
    GET_INPUT_VARIANTS(exposure_input_data)

    // run MR and coloc for both cis and all
    CIS(GET_INPUT_VARIANTS.out.cis_exposure, GET_INPUT_VARIANTS.out.cis_outcome, GET_INPUT_VARIANTS.out.cis_harm)
    ALL(GET_INPUT_VARIANTS.out.top_exposure, GET_INPUT_VARIANTS.out.top_outcome, GET_INPUT_VARIANTS.out.all_harm)
}
