
// import local modules
include { GET_INPUT_VARIANTS } from '../modules/get_input_variants.nf'

// main pipeline
workflow MULTI_EXPOSURE_MR {

    // get exposure file names
    exposure_input_data = Channel.fromPath("${params.exposure_input_dir}/*${params.data_format.exposure}")

    // reformat data into TwoSampleMR format w/ significant variants
    GET_INPUT_VARIANTS( exposure_input_data )

    // GET_INPUT_VARIANTS.out.tsmr_outcome.view()

    // do plink clumping
    // INSTRUMENT_SELECTION

    // get ld matrix
    // EXTRACT_LD_MATRIX

    // run mr
    // RUN_MR

    // do coloc for cis region and top mr variant (plus plot these regions)
    // RUN_COLOC

    // generate mr plots
    // PLOT_MR
}
