
// import local modules
include { GET_INPUT_VARIANTS } from '../modules/get_input_variants.nf'
include { GET_LD_MATRIX } from '../modules/get_ld_matrix.nf'

// main pipeline
workflow MULTI_EXPOSURE_MR {

    // get exposure file names
    exposure_input_data = Channel.fromPath("${params.exposure.input_dir}/*")

    // reformat data into TwoSampleMR format w/ significant variants & clump
    GET_INPUT_VARIANTS(exposure_input_data)

    // do coloc for cis region and top mr variant (plus plot these regions)
    // RUN_COLOC(exposure_input_data, GET_INPUT_VARIANTS.out.cis_location)
    // RUN_COLOC(exposure_input_data, GET_INPUT_VARIANTS.out.top_location)

    // run mrs
    // MR_WORKFLOW(GET_INPUT_VARIANTS.out.all_harm)
    // MR_WORKFLOW(GET_INPUT_VARIANTS.out.cis_harm)
}
