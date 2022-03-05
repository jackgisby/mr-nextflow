
// import local modules
include { RUN_COLOC } from "$baseDir/modules/run_coloc.nf"
include { RUN_MR } from "$baseDir/modules/run_mr.nf"
include { COLLATE_MR } from "$baseDir/modules/collate_results.nf"
include { COLLATE_COLOC } from "$baseDir/modules/collate_results.nf"

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

    COLLATE_MR(
        RUN_MR.out.mr_results.collectFile(name: "mr_results.csv", keepHeader: true, skip: 1, newLine: false),
        RUN_MR.out.singlesnp.collectFile(name: "singlesnp.csv", keepHeader: true, skip: 1, newLine: false),
        RUN_MR.out.directionality.collectFile(name: "directionality.csv", keepHeader: true, skip: 1, newLine: false),
        RUN_MR.out.heterogeneity.collectFile(name: "heterogeneity.csv", keepHeader: true, skip: 1, newLine: false),
        RUN_MR.out.pleiotropy.collectFile(name: "pleiotropy.csv", keepHeader: true, skip: 1, newLine: false),
        RUN_MR.out.mr_results_leaveoneout.collectFile(name: "mr_results_leaveoneout.csv", keepHeader: true, skip: 1, newLine: false)
    )

    if (params.coloc.run) {
        RUN_COLOC(exposure, outcome, LD)

        COLLATE_COLOC(
            RUN_COLOC.out.coloc_res.collectFile(name: "coloc_res.csv", keepHeader: true, skip: 1, newLine: false),
            RUN_COLOC.out.coloc_susie_res.collectFile(name: "coloc_susie_res.csv", keepHeader: true, skip: 1, newLine: false)
        )
    }
}
