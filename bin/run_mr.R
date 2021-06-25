# load libraries
library("optparse")
library("TwoSampleMR")

# do argument parsing
option_list = list(
    make_option(c("--harm"), type="character", default=NULL, help="harmonised dataset in TwoSampleMR format", metavar="character"),
    make_option(c("--mrinput"), type="character", default=NULL, help="harmonised dataset in MendelianRandomization format", metavar="character"),
    make_option(c("--auxiliary_script_dir"), type="character", default=NULL, help="the location of helper scripts", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

print("options")
print(opt)

source(paste0(opt$auxiliary_script_dir, "/mr.R"))

get_blank_return <- function() {

    result_colnames <- list(
        "mr_results"=c("id.exposure", "id.outcome", "outcome", "exposure", "method", "nsnp", "b", "se", "pval", "package"), 
        "singlesnp"=c("exposure", "outcome", "id.exposure", "id.outcome", "samplesize", "SNP", "b", "se", "p"), 
        "directionality"=c("id.exposure", "id.outcome", "exposure", "outcome", "snp_r2.exposure", "snp_r2.outcome", "correct_causal_direction", "steiger_pval", "package"), 
        "heterogeneity"=c("id.exposure", "id.outcome", "outcome", "exposure", "method", "Q", "Q_df", "Q_pval", "package"), 
        "pleiotropy"=c("id.exposure", "id.outcome", "outcome", "exposure", "egger_intercept", "egger_se", "egger_pval", "package")
    )

    blank_results <- list()

    for (result_type in names(result_colnames)) {
        
        blank_results[[result_type]] <- data.frame(matrix(ncol = length(result_colnames[[result_type]]), nrow = 0))
        colnames(blank_results[[result_type]]) <- result_colnames[[result_type]]
    }

    return(blank_results)
}

run_mr <- function(opt) {

    exposure_name <- gsub(".csv", "", gsub("_harm", "", opt$harm))

    harm <- data.frame(readr::read_csv(opt$harm))
    print("harm")
    print(head(harm))

    mrinput <- readRDS(opt$mrinput)
    print("mrinput")
    print(names(mrinput))

    return_names <- c(paste0(exposure_name, "_mr_results"), paste0(exposure_name, "_singlesnp"), paste0(exposure_name, "_directionality"), paste0(exposure_name, "_heterogeneity"), paste0(exposure_name, "_pleiotropy"))
    null_return <- get_blank_return()
    names(null_return) <- return_names

    if (nrow(harm) == 0) {
        return(null_return)
    }

    mr_parameters <- default_parameters()
    method_list <- c('mr_wald_ratio', 'mr_egger_regression', 'mr_ivw', 'mr_ivw_fe', 'mr_weighted_median', 'mr_weighted_mode', 'mr_two_sample_ml')

    # generate MR results with TwoSampleMR
    mr_results <- mr(harm, parameters = mr_parameters, method_list = method_list)
    mr_results$package <- "TwoSampleMR"
    
    # add MendelianRandomization MR results and sensitivity tests
    mr_results <- add_tests(mr_results, harm, mrinput, mr_parameters, method_list)

    # set up final results
    names(mr_results) <- return_names
    return(mr_results)
}

files_out <- run_mr(opt)

for (i in 1:length(files_out)) {
    write.csv(files_out[[i]], paste0(names(files_out)[i], ".csv"), row.names = FALSE)
}
