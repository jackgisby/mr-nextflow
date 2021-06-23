
// import local modules
include { RUN_MR } from '../modules/run_mr.nf'
include { COLLATE_RESULTS } from '../modules/collate_results.nf'

// main pipeline
workflow MR_WORKFLOW {

    take: harm

    main:
    // get ld matrix if required
    if (params.use_ld_matrix) {

        GET_LD_MATRIX(harm)

        GET_LD_MATRIX.out.mr_input.to{mr_input}
        GET_LD_MATRIX.out.matrix.to{matrix}
        
    } else {
        Channel.value("").to{mr_input}
        Channel.value("").to{matrix}
    }

    RUN_MR(harm, mr_input, matrix)
    COLLATE_RESULTS(RUN_MR.out.mr_results)
}
