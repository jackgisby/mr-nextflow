
// import local modules
include { RUN_COLOC } from "$baseDir/modules/run_coloc.nf"
include { RUN_MR } from "$baseDir/modules/run_mr.nf"

// run mr and coloc - do for both cis and all
workflow CIS {

    take:
    exposure
    outcome
    LD
    harm
    mrinput

    main:
    MR_COLOC(exposure, outcome, LD, harm, mrinput)
}

workflow ALL {

    take:
    exposure
    outcome
    LD
    harm
    mrinput

    main:  // same workflow, but can't run same thing twice with same name
    MR_COLOC(exposure, outcome, LD, harm, mrinput)
}

workflow MR_COLOC {

    take:
    exposure
    outcome
    LD
    harm
    mrinput

    main:
    RUN_MR(harm, mrinput)
    RUN_COLOC(exposure, outcome, LD)
    // COLLATE_MR
    // COLLATE_COLOC
}
