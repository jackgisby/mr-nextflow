
// import local modules
include { RUN_COLOC } from "$baseDir/modules/run_coloc.nf"
include { RUN_MR } from "$baseDir/modules/run_mr.nf"

// run mr and coloc - do for both cis and all
workflow CIS {

    take:
    exposure
    outcome
    harm

    main:
    MR_COLOC(exposure, outcome, harm)
}

workflow ALL {

    take:
    exposure
    outcome
    harm

    main:  // same workflow, but can't run same thing twice with same name
    MR_COLOC(exposure, outcome, harm)
}

workflow MR_COLOC {

    take:
    exposure
    outcome
    harm

    main:
    // RUN_MR(harm)
    RUN_COLOC(exposure, outcome)
}
