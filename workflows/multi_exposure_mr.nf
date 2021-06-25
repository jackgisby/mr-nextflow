// import local modules
include { GET_INPUT_VARIANTS } from "$baseDir/modules/get_input_variants.nf"
include { CIS } from "$baseDir/workflows/workflow_split.nf"
include { ALL } from "$baseDir/workflows/workflow_split.nf"

// main pipeline
workflow MULTI_EXPOSURE_MR {

    // get exposure file names
    exposure_input_data = Channel.fromPath("${params.exposure.input_dir}/*")

    // reformat data into TwoSampleMR format w/ significant variants & clump
    GET_INPUT_VARIANTS(exposure_input_data)  // TODO: collectFile or similar to remove empty runs?

    // run MR and coloc for both cis and all - these workflows do the same thing
    CIS(GET_INPUT_VARIANTS.out.cis_exposure, GET_INPUT_VARIANTS.out.cis_outcome, GET_INPUT_VARIANTS.out.cis_LD, GET_INPUT_VARIANTS.out.cis_harm, GET_INPUT_VARIANTS.out.cis_mrinput)
    ALL(GET_INPUT_VARIANTS.out.top_exposure, GET_INPUT_VARIANTS.out.top_outcome, GET_INPUT_VARIANTS.out.top_LD, GET_INPUT_VARIANTS.out.all_harm, GET_INPUT_VARIANTS.out.all_mrinput)
}
