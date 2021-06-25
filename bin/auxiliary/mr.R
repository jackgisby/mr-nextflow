add_tests <- function(mr_results, harm, mrinput, mr_parameters, method_list) {
    
    if (nrow(mr_results) <= 1) {
        return(mr_results)
    }
    
    # mr ivw random
    ivw_cor <- MendelianRandomization::mr_ivw(
        mrinput,
        model = "random",
        robust = FALSE,
        penalized = FALSE,
        correl = TRUE,
        weights = "simple"
    )
    
    mr_results <- rbind(mr_results, data.frame(
        id.exposure = harm$id.exposure[1], 
        id.outcome = harm$id.outcome[1], 
        outcome = mrinput@outcome, 
        exposure = mrinput@exposure, 
        method = "Inverse variance weighted", 
        nsnp = length(mrinput@snps), 
        b = ivw_cor@Estimate, 
        se = ivw_cor@StdError, 
        pval = ivw_cor@Pvalue,
        package = "MendelianRandomization"
    ))
    
    # mr ivw fixed
    ivw_fe_cor <- MendelianRandomization::mr_ivw(
        mrinput,
        model = "fixed",
        robust = FALSE,
        penalized = FALSE,
        correl = TRUE,
        weights = "simple"
    )
    
    mr_results <- rbind(mr_results, data.frame(
        id.exposure = harm$id.exposure[1], 
        id.outcome = harm$id.outcome[1], 
        outcome = mrinput@outcome, 
        exposure = mrinput@exposure, 
        method = "Inverse variance weighted (fixed effects)", 
        nsnp = length(mrinput@snps), 
        b = ivw_fe_cor@Estimate, 
        se = ivw_fe_cor@StdError, 
        pval = ivw_fe_cor@Pvalue,
        package = "MendelianRandomization"
    ))
    
    # mr_maxlik
    ml_cor <- MendelianRandomization::mr_maxlik(
        mrinput, 
        correl = TRUE
    )
    
    mr_results <- rbind(mr_results, data.frame(
        id.exposure = harm$id.exposure[1], 
        id.outcome = harm$id.outcome[1], 
        outcome = mrinput@outcome, 
        exposure = mrinput@exposure, 
        method = "Maximum likelihood", 
        nsnp = length(mrinput@snps), 
        b = ml_cor@Estimate,
        se = ml_cor@StdError, 
        pval = ml_cor@Pvalue,
        package = "MendelianRandomization"
    ))

    heterogeneity <- mr_heterogeneity(
        harm, 
        parameters = mr_raps_parameters
    )
    heterogeneity$package <- "TwoSampleMR"
    
    heterogeneity_cor <- data.frame(
        id.exposure = harm$id.exposure[1],
        id.outcome = harm$id.outcome[1],
        outcome = mrinput@outcome,
        exposure = mrinput@exposure,
        method = "Inverse variance weighted",
        Q = ivw_cor@Heter.Stat[1],
        Q_df = NA,
        Q_pval = ivw_cor@Heter.Stat[2],
        package = "MendelianRandomization"
    )

    heterogeneity <- rbind(heterogeneity, heterogeneity_cor)

    singlesnp <- mr_singlesnp(
        harm, 
        parameters = mr_parameters,
        single_method = "mr_wald_ratio",
        all_method = c("mr_ivw")
    )
    
    directionality_test <- directionality_test(harm)  # only for tsmr package
    directionality_test$package <- "TwoSampleMR"

    pleiotropy_cols <- c("id.exposure", "id.outcome", "outcome", "exposure", "egger_intercept", "egger_se", "egger_pval", "package")
    pleiotropy <- data.frame(matrix(ncol = length(pleiotropy_cols), nrow=0))  # only run if >2 SNPs but need to have right colnames either way
    colnames(pleiotropy) <- pleiotropy_cols

    
    if (length(mrinput@snps) > 2) {
        
        # egger
        egger_cor <- MendelianRandomization::mr_egger(
            mrinput, 
            robust = FALSE,
            penalized = FALSE,
            correl = TRUE
        )
        
        mr_results <- rbind(mr_results, data.frame(
            id.exposure = harm$id.exposure[1], 
            id.outcome = harm$id.outcome[1], 
            outcome = mrinput@outcome, 
            exposure = mrinput@exposure, 
            method = "MR Egger", 
            nsnp = length(mrinput@snps), 
            b = egger_cor@Estimate,
            se = egger_cor@StdError.Est, 
            pval = egger_cor@Pvalue.Est,
            package = "MendelianRandomization"
        ))

        pleiotropy <-  mr_pleiotropy_test(harm)
        pleiotropy$package="TwoSampleMR"

        pleiotropy_cor <- data.frame(
            id.exposure = harm$id.exposure[1],
            id.outcome = harm$id.outcome[1],
            outcome = mrinput@outcome,
            exposure = mrinput@exposure,
            egger_intercept = egger_cor@Intercept,
            se = egger_cor@StdError.Int,
            pval = egger_cor@Pleio.pval,
            package="MendelianRandomization"
        )

        pleiotropy <- rbind(pleiotropy, pleiotropy_cor)
        colnames(pleiotropy) <- sapply(colnames(pleiotropy), function(cname) {
            if (cname == "se") {
                return("egger_se")
            } else if (cname == "pval") {
                return("egger_pval")
            } else {
                return(cname)
            }
        })
        
        # mr_median
        weighted_median_cor <- MendelianRandomization::mr_median(
            mrinput, 
            weighting = "weighted"
        )
        
        mr_results <- rbind(mr_results, data.frame(
            id.exposure = harm$id.exposure[1], 
            id.outcome = harm$id.outcome[1], 
            outcome = mrinput@outcome, 
            exposure = mrinput@exposure, 
            method = "Weighted median", 
            nsnp = length(mrinput@snps), 
            b = weighted_median_cor@Estimate,
            se = weighted_median_cor@StdError, 
            pval = weighted_median_cor@Pvalue,
            package = "MendelianRandomization"
        ))
        
        # mr_mbe
        weighted_mode_cor <- MendelianRandomization::mr_median(
            mrinput, 
            weighting = "weighted"
        )
        
        mr_results <- rbind(mr_results, data.frame(
            id.exposure = harm$id.exposure[1], 
            id.outcome = harm$id.outcome[1], 
            outcome = mrinput@outcome, 
            exposure = mrinput@exposure, 
            method = "Weighted mode", 
            nsnp = length(mrinput@snps), 
            b = weighted_mode_cor@Estimate,
            se = weighted_mode_cor@StdError, 
            pval = weighted_mode_cor@Pvalue,
            package = "MendelianRandomization"
        ))
    }

    return(list(
        "mr_results"=mr_results, 
        "singlesnp"=singlesnp, 
        "directionality"=directionality_test,
        "heterogeneity"=heterogeneity,
        "pleiotropy"=pleiotropy
    ))
}

get_collated_df <- function(mr_results, num_exposures, directionality, heterogeneity, pleiotropy) {

    collated <- mr_results
    sensitivity_join_columns <- c("id.exposure"="id.exposure", "id.outcome"="id.outcome", "outcome"="outcome", "exposure"="exposure", "package"="package")

    collated <- dplyr::left_join(collated, directionality, by=sensitivity_join_columns)
    collated <- dplyr::left_join(collated, heterogeneity, by=c(sensitivity_join_columns, "method"))

    if (nrow(pleiotropy) > 0) {
        collated <- dplyr::left_join(collated, pleiotropy, by=sensitivity_join_columns)
    }

    stopifnot(nrow(mr_results) == nrow(collated))
    
    collated <- collated[collated$method %in% c("Inverse variance weighted", "Wald ratio"),]
    collated <- collated[!(collated$nsnp == 1 & collated$method != "Wald ratio"),]

    collated$pval_bh <- 1
    collated$pval_bonf <- 1

    for (method in unique(mr_results$method)) {
        collated <- add_method(collated, mr_results, method, sensitivity_join_columns)
    }

    for (package in unique(collated$package)) {

        collated$pval_bh[collated$package == package] <- p.adjust(collated$pval[collated$package == package], method = "BH")
        collated$pval_bonf[collated$package == package] <- p.adjust(collated$pval[collated$package == package], method = "bonferroni", n = num_exposures)

        exposures <- collated$exposure[collated$package == package]
        stopifnot(length(exposures) == length(unique(exposures)))
    }

    return(collated)
}

# add method as a column
add_method <- function(collated, full_mr_results, method, join_columns) {

    if (!(method %in% full_mr_results$method)) {
        return(collated)
    }
    
    results_to_add <- full_mr_results[full_mr_results$method == method,]
    results_to_add <- dplyr::select(results_to_add, "id.exposure", "id.outcome", "exposure", "outcome", "package", "b", "se", "pval")
    colnames(results_to_add) <- c("id.exposure", "id.outcome", "exposure", "outcome", "package", paste0("b_", method), paste0("se_", method), paste0("pval_", method))

    collated <- dplyr::left_join(collated, results_to_add, by=join_columns)
    return(collated)
}
