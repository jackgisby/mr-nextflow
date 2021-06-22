#' Function for applying to MR to a given exposure and outcome
#' 
#' @param exposure_name
#' Name of the exposure - e.g. gene name corresponding to a QTL. Used to find
#' the gene position in `gene_position_map_file` if this is provided.
#' 
#' @param exposure_gwas
#' A `data.frame` in the format returned by \code{\link[TwoSampleMR]{format_data}}.
#' Contains the variants for the exposure, for instance a QTL.
#' 
#' @param outcome_gwas
#' A `data.frame` in the format returned by \code{\link[TwoSampleMR]{format_data}}.
#' Contains the variants for the outcome.
#' 
#' @param linkage_file
#' File(s) for the calculation of linkage disequilibrium. If this is a single file,
#' this should be a string referring to the path of the file. If there is one file
#' per chromosome, then this should be a list; the keys of the list should be
#' integers referring to the chromosome number and the values should be a string
#' referring to the path of the file.
#' 
#' @param p_threshold 
#' Maxmimum p value threshold for variants in `exposure_gwas`.
#' 
#' @param gene_position_map_file
#' If `NULL`, only trans exposure variants are tested. If given, should be a 
#' TSV containing the following columns:
#' \itemize{
#'     \item chromosome_name
#'     \item start_position
#'     \item end_position
#'     \item ensembl_gene_id
#'     \item hgnc_symbol
#' }
#' 
#' `exposure_name` will be used to look up the `?` column and the corresponding
#' locus will be used to define which exposure variants are cis.
#' 
#' @param clump_kb
#' Parameter clump-kb for PLINK clumping. 
#' 
#' @param clump_r2
#' Parameter clump-r2 for PLINK clumping. 
#' 
#' @param cis_region
#' The threshold around a gene for which to consider exposure variants as cis.
#' 
#' @param sink_plink
#' If TRUE, `sink(/dev/null)` will be called around the call to PLINK.
#' 
#' @param linkage_file_format
#' The file format to be used as input to PLINK.
#' 
#' @return 
#' A large list of results.
#' 
#' @export

mr_workflow <- function(
    exposure_name,
    exposure_gwas,
    outcome_gwas,
    linkage_file,
    p_threshold = 5e-8,
    results_dir = NULL,
    gene_position_map_file = NULL,
    clump_kb = 10000,
    clump_r2 = 0.01,
    cis_region = 1000000,
    sink_plink = TRUE,
    plink_bin = "plink",
    wald_only = FALSE
) {
    if (!is.null(results_dir)) {
        
        dir.create(results_dir, showWarnings = FALSE)
        
        if (length(list.files(results_dir)) > 0) {
            message("Warning: results directory is not empty; previous results may be overwritten or mixed with current results.")
        }
        
        plot_dir <- file.path(results_dir, "plots")
        
    } else {
        plot_dir <- NULL
    }

    mr_results <- list()
    mr_results[["run_info"]] <- list()
    
    mr_results[["run_info"]][["exposure_name"]] <- exposure_name
    mr_results[["run_info"]][["linkage_file"]] <- linkage_file
    mr_results[["run_info"]][["p_threshold"]] <- p_threshold
    mr_results[["run_info"]][["clump_kb"]] <- clump_kb
    mr_results[["run_info"]][["clump_r2"]] <- clump_r2
    mr_results[["run_info"]][["cis_region"]] <- cis_region
    
    #----------------- prune LD
    
    if (any(exposure_gwas$pval.exposure <= p_threshold)) {
        
        exposure_gwas <- exposure_gwas[which(exposure_gwas$pval.exposure <= p_threshold),]
        
    } else {
        print("No significant exposure SNPs")
        return()
    }
    
    if (nrow(exposure_gwas) > 1) {
        
        clump <- plink_clump(exposure_gwas, clump_kb = clump_kb, clump_r2 = clump_r2, clump_p = 1,
                             linkage_file = linkage_file, sink_plink = sink_plink,
                             plink_bin = plink_bin)
        
        before_clump <- nrow(exposure_gwas)
        exposure_gwas <- exposure_gwas[which(exposure_gwas$SNP %in% clump$SNP),]
        
        if (nrow(exposure_gwas) == 0) {
            
            print("No significant exposure SNPs found following clumping")
            print(paste0(before_clump, " SNPs found prior to clumping."))
            return()
        }
    } 
    
    exposure_gwas$type <- "trans"
    #----------------- define genes/loci
    
    # get gene position to define cis loci
    if (!is.null(gene_position_map_file)) {
        
        gene_position_map <- fread(gene_position_map_file)
        gene_position_map <- gene_position_map[which(gene_position_map$chromosome_name %in% 1:22),]
        act_map <- gene_position_map[tolower(gene_position_map$hgnc_symbol) == exposure_name,]
        
        # define which locus is a cis-pqtl if any - variant within gene/cis region
        if (nrow(act_map) > 0) {
            
            exposure_gwas$type[exposure_gwas$chr.exposure == act_map$chromosome_name &
                                   exposure_gwas$pos.exposure > act_map$start_position - cis_region & 
                                   exposure_gwas$pos.exposure < act_map$end_position + cis_region
            ] <- "cis"
        }
        
        mr_results[["run_info"]][["inc_cis_trans"]] <- TRUE
        
    } else {
        message("gene_position_map_file is NULL; will find trans exposure variants only")
        
        mr_results[["run_info"]][["inc_cis_trans"]] <- FALSE
    }
    
    if (any(is.na(exposure_gwas$type))) {
        
        message("NA type for a clumped variant")
        return()
    }
    
    #----------------- wald only
    
    if (wald_only) {
        exposure_gwas <- exposure_gwas[exposure_gwas$type == "cis",]
        exposure_gwas <- exposure_gwas[which.min(exposure_gwas$pval.exposure),]
    }
    
    if (nrow(exposure_gwas) == 0) {
        
        print("No significant instruments found after restriction to wald only")
        return()
    }
    
    #----------------- remove MHC
    
    nrow_prior_to_mhc_removal <- nrow(exposure_gwas)
    MHC <- c(26000000, 34000000)
    mr_results[["run_info"]][["mhc"]] <- MHC
    
    bad_mhc <- which(exposure_gwas$chr.exposure == 6 &
                         exposure_gwas$pos.exposure >= MHC[1] &
                         exposure_gwas$pos.exposure <= MHC[2])
    
    if (length(bad_mhc) > 0) {
        exposure_gwas <- exposure_gwas[-bad_mhc,]
    }
    
    if (nrow(exposure_gwas) == 0) {
        
        print("No significant instruments found after removal of MHC region")
        print(paste0(nrow_prior_to_mhc_removal, " significant SNPs were removed"))
        return()
    }
    
    #----------------- run mrs
    
    mr_parameters <- default_parameters()
    method_list <- c('mr_wald_ratio', 'mr_egger_regression', 'mr_ivw', 'mr_ivw_fe', 'mr_weighted_median', 'mr_weighted_mode', 'mr_two_sample_ml')
    
    mr_results[["run_info"]][["method_list"]] <- method_list
    mr_results[["run_info"]][["which_mr_runs"]] <- c()
    mr_results[["run_info"]][["mr_parameters"]] <- mr_parameters
    
    # cis data
    if (any(exposure_gwas$type == "cis")) {
        
        mr_results[["cis_harm"]] <- exposure_gwas[which(exposure_gwas$type == "cis"),]
        mr_results[["cis_harm"]]$exposure <- paste(exposure_name, "cis", sep = "_")
        mr_results[["cis_harm"]] <- harmonise_data(mr_results[["cis_harm"]], outcome_gwas, 1)
        
        mr_results[["run_info"]][["cis_num_removed_outliers"]] <- 0
        if (nrow(mr_results[["cis_harm"]]) > 2) {
            mr_results[["cis_harm_removed_outliers"]] <- try(remove_outliers_radial(mr_results[["cis_harm"]]))
            
            if (inherits(mr_results[["cis_harm_removed_outliers"]], "data.frame")) {
                mr_results[["run_info"]][["cis_num_removed_outliers"]] <- nrow(mr_results[["cis_harm"]]) - nrow(mr_results[["cis_harm_removed_outliers"]])
            }        
        }
        
        mr_results[["cis_mr"]] <- mr(mr_results[["cis_harm"]], 
                                     parameters = mr_results[["run_info"]][["mr_parameters"]], 
                                     method_list = mr_results[["run_info"]][["method_list"]])
        mr_results[["cis_mr"]]$package <- "TwoSampleMR"
        
        mr_results <- ld_cor_mr(mr_results, exposure_type = "cis", linkage_file, sink_plink, plink_bin)
        
        mr_results[["run_info"]][["which_mr_runs"]] <- c(mr_results[["run_info"]][["which_mr_runs"]], "cis") 
    }
    
    # trans data 
    if (any(exposure_gwas$type == "trans")) {
        mr_results[["trans_harm"]] <- exposure_gwas[which(exposure_gwas$type == "trans"),]
        mr_results[["trans_harm"]]$exposure <- paste(exposure_name, "trans", sep = "_")
        mr_results[["trans_harm"]] <- harmonise_data(mr_results[["trans_harm"]], outcome_gwas, 1)
        
        mr_results[["run_info"]][["trans_num_removed_outliers"]] <- 0
        if (nrow(mr_results[["trans_harm"]]) > 2) {
            mr_results[["trans_harm_removed_outliers"]] <- try(remove_outliers_radial(mr_results[["trans_harm"]]))
            
            if (inherits(mr_results[["trans_harm_removed_outliers"]], "data.frame")) {
                mr_results[["run_info"]][["trans_num_removed_outliers"]] <- nrow(mr_results[["trans_harm"]]) - nrow(mr_results[["trans_harm_removed_outliers"]])
            }        
        }
        
        mr_results[["trans_mr"]] <- mr(mr_results[["trans_harm"]], 
                                       parameters = mr_results[["run_info"]][["mr_parameters"]], 
                                       method_list = mr_results[["run_info"]][["method_list"]])
        mr_results[["trans_mr"]]$package <- "TwoSampleMR"
        
        mr_results <- ld_cor_mr(mr_results, exposure_type = "trans", linkage_file, sink_plink, plink_bin)
        
        mr_results[["run_info"]][["which_mr_runs"]] <- c(mr_results[["run_info"]][["which_mr_runs"]], "trans") 
    }

    # all data
    exposure_gwas$exposure <- paste(exposure_name, "all", sep = "_")
    mr_results[["all_harm"]] <- harmonise_data(exposure_gwas, outcome_gwas, 1)
    
    mr_results[["run_info"]][["all_num_removed_outliers"]] <- 0
    if (nrow(mr_results[["all_harm"]]) > 2) {
        mr_results[["all_harm_removed_outliers"]] <- try(remove_outliers_radial(mr_results[["all_harm"]]))

        if (inherits(mr_results[["all_harm_removed_outliers"]], "data.frame")) {
            mr_results[["run_info"]][["all_num_removed_outliers"]] <- nrow(mr_results[["all_harm"]]) - nrow(mr_results[["all_harm_removed_outliers"]])
        }
    }
    
    mr_results[["run_info"]][["which_mr_runs"]] <- c(mr_results[["run_info"]][["which_mr_runs"]], "all") 
    mr_results[["all_mr"]] <- mr(mr_results[["all_harm"]], 
                                 parameters = mr_results[["run_info"]][["mr_parameters"]], 
                                 method_list = mr_results[["run_info"]][["method_list"]])
    mr_results[["all_mr"]]$package <- "TwoSampleMR"
    
    mr_results <- ld_cor_mr(mr_results, exposure_type = "all", linkage_file, sink_plink, plink_bin)
    
    # extra analyses
    mr_results[["run_info"]][["additional_analyses"]] <- c('mr', 'harm', 'mr_heterogeneity', 'mr_singlesnp', 'mr_leaveoneout', 'pleiotropy_test', 'directionality_test', 'mr_raps', 'ivw_radial')
    mr_results <- run_additional_analyses(mr_results)
    
    # plots
    if (!is.null(plot_dir) & nrow(exposure_gwas) > 1) {
        
        dir.create(results_dir, showWarnings = FALSE)
        
        plot_dir <- file.path(results_dir, "plots")
        dir.create(plot_dir, showWarnings = FALSE)
        
        
        mr_results[["run_info"]][["plots"]] <- TRUE
        mr_results[["plots"]] <- make_mr_plots(mr_results, plot_dir)
    }
    
    #----------------- save files
    
    if (!is.null(results_dir)) {
        
        # save mr_results R object and main results as tsv
        saveRDS(mr_results, paste0(results_dir, "/results.rds"))
        
        collated_results_df <- c()
        
        # add to results
        for (which_mr_run in mr_results[["run_info"]][["which_mr_runs"]]) {
            
            collated_results_df <- rbind(collated_results_df, mr_results[[paste0(which_mr_run, "_mr")]])
        }
        
        write.table(collated_results_df, 
                    file = paste0(results_dir, "_mr_results.tsv"), 
                    row.names = F, quote = F, sep = "\t")
        
        # save run details
        sink(paste0(results_dir, "/run_info.txt"))
        print(mr_results[["run_info"]])
        sink()
    }
    
    # return results for collation and downstream saving
    return(mr_results)
}

#' runs additional analyses based on the main pipeline results
#' 
#' For each mr run (i.e. for each of cis, trans, all) performs additional
#' analyses as specified in the workflow.
#' 
#' This function is just an extension of `mr_workflow`. The analyses that
#' are to be carried out here are specified within 
#' `mr_results[["run_info"]][["additional_analyses"]]`.
#' 
#' @param mr_results
#' A large list generated within `mr_workflow`
#' 
#' @return 
#' Returns the original list, with additional entries.

run_additional_analyses <- function(mr_results) {

    for (which_mr_run in mr_results[["run_info"]][["which_mr_runs"]]) {
        
        if ("mr_heterogeneity" %in% mr_results[["run_info"]][["additional_analyses"]] & nrow(mr_results[[paste0(which_mr_run, "_harm")]]) >= 2) {
            
            h_method <- subset(mr_method_list(), heterogeneity_test)$obj
            h_method <- h_method[h_method %in% mr_results[["run_info"]]["method_list"][[1]]]
            
            mr_results[[paste0(which_mr_run, "_mr_heterogeneity")]] <- mr_heterogeneity(
                mr_results[[paste0(which_mr_run, "_harm")]], 
                parameters = mr_results[["run_info"]][["mr_parameters"]],
                method_list = h_method
            )
            
            mr_results[[paste0(which_mr_run, "_mr_heterogeneity_cor")]] <- data.frame(
                id.exposure = mr_results[[paste0(which_mr_run, "_harm")]]$id.exposure[1],
                id.outcome = mr_results[[paste0(which_mr_run, "_harm")]]$id.outcome[1],
                outcome = mr_results[[paste0(which_mr_run, "_mr_object")]]@outcome,
                exposure = mr_results[[paste0(which_mr_run, "_mr_object")]]@exposure,
                method = "Inverse variance weighted",
                Q = mr_results[[paste0(which_mr_run, "_ivw_cor")]]@Heter.Stat[1],
                Q_df = NA,
                Q_pval = mr_results[[paste0(which_mr_run, "_ivw_cor")]]@Heter.Stat[2]
            )
        }
        
        if ("pleiotropy_test" %in% mr_results[["run_info"]][["additional_analyses"]] & nrow(mr_results[[paste0(which_mr_run, "_harm")]]) >= 2) {
            
            mr_results[[paste0(which_mr_run, "_pleiotropy_test")]] <-  mr_pleiotropy_test(
                mr_results[[paste0(which_mr_run, "_harm")]]
            )
            
            if (!is.null(mr_results[[paste0(which_mr_run, "_mr_egger_regression_cor")]])) {
                mr_results[[paste0(which_mr_run, "_pleiotropy_test_cor")]] <- data.frame(
                    id.exposure = mr_results[[paste0(which_mr_run, "_harm")]]$id.exposure[1],
                    id.outcome = mr_results[[paste0(which_mr_run, "_harm")]]$id.outcome[1],
                    outcome = mr_results[[paste0(which_mr_run, "_mr_object")]]@outcome,
                    exposure = mr_results[[paste0(which_mr_run, "_mr_object")]]@exposure,
                    egger_intercept = mr_results[[paste0(which_mr_run, "_mr_egger_regression_cor")]]@Intercept,
                    se = mr_results[[paste0(which_mr_run, "_mr_egger_regression_cor")]]@StdError.Int,
                    pval = mr_results[[paste0(which_mr_run, "_mr_egger_regression_cor")]]@Pleio.pval
                )
            }
        }
        
        if ("mr_singlesnp" %in% mr_results[["run_info"]][["additional_analyses"]]) {
            
            mr_results[[paste0(which_mr_run, "_mr_singlesnp")]] <- mr_singlesnp(
                mr_results[[paste0(which_mr_run, "_harm")]], 
                parameters = mr_results[["run_info"]][["mr_parameters"]],
                single_method = "mr_wald_ratio",
                all_method = mr_results[["run_info"]][["method_list"]]
            )
        }
        
        if ("mr_leaveoneout" %in% mr_results[["run_info"]][["additional_analyses"]]) {
            
            mr_results[[paste0(which_mr_run, "_mr_leaveoneout")]] <- mr_leaveoneout(
                mr_results[[paste0(which_mr_run, "_harm")]], 
                parameters = mr_results[["run_info"]][["mr_parameters"]]
            )
        }
        
        if ("directionality_test" %in% mr_results[["run_info"]][["additional_analyses"]]) {
            
            mr_results[[paste0(which_mr_run, "_directionality_test")]] <- directionality_test(
                mr_results[[paste0(which_mr_run, "_harm")]]
            )
        }
        
        if ("ivw_radial" %in% mr_results[["run_info"]][["additional_analyses"]] & nrow(mr_results[[paste0(which_mr_run, "_harm")]]) >= 2) {

            mr_results[[paste0(which_mr_run, "_ivw_radial")]] <- get_ivw_radial_results(
                mr_results[[paste0(which_mr_run, "_harm")]]
            )
        }
        
        # run raps with different parameters
        mr_results[["run_info"]][["mr_raps_parameters"]] <- mr_results[["run_info"]][["mr_parameters"]]
        mr_results[["run_info"]][["mr_raps_parameters"]][["over.dispersion"]] <- FALSE
        
        mr_raps_results <- mr(
            mr_results[[paste0(which_mr_run, "_harm")]], 
            parameters = mr_results[["run_info"]][["mr_raps_parameters"]],
            method_list = "mr_raps"
        )
        mr_raps_results$package <- "TwoSampleMR"
        
        mr_results[[paste0(which_mr_run, "_mr")]] <- rbind(mr_results[[paste0(which_mr_run, "_mr")]], mr_raps_results)
    }
    
    return(mr_results)
}

#' Makes plots for the MR workflow
#' 
#' Takes the large MR result list from `mr_workflow` and generates plots.
#' 
#' @param mr_results
#' The large list from `mr_workflow`.
#' 
#' @param plot_dir
#' The directory in which to save plots to. If NULL, will not save plots.
#' 
#' @param print_plots
#' If TRUE, prints the generated plots in addition to saving them.

make_mr_plots <- function(
    mr_results,             
    plot_dir=NULL,
    print_plots=FALSE
) {
    plots <- list()
    
    for (which_mr_run in mr_results[["run_info"]][["which_mr_runs"]]) {
        
        plots[[paste0(which_mr_run, "_scatter")]] <- mr_scatter_plot(mr_results[[paste0(which_mr_run, "_mr")]], mr_results[[paste0(which_mr_run, "_harm")]])
        
        plots[[paste0(which_mr_run, "_scatter")]][[1]] <- plots[[paste0(which_mr_run, "_scatter")]][[1]] +
            geom_hline(yintercept = 0) +
            geom_vline(xintercept = 0)
        
        if (length(mr_results[["run_info"]][["which_mr_runs"]]) > 1 & which_mr_run == "all") {
            
            plots[[paste0(which_mr_run, "_scatter")]][[1]] <- plots[[paste0(which_mr_run, "_scatter")]][[1]] + 
                geom_point(aes(col = type))
        }
        
        if ("mr_singlesnp" %in% mr_results[["run_info"]][["additional_analyses"]]) {
            plots[[paste0(which_mr_run, "_mr_singlesnp_forest")]] <- mr_forest_plot(mr_results[[paste0(which_mr_run, "_mr_singlesnp")]])
            plots[[paste0(which_mr_run, "_mr_singlesnp_funnel")]] <- mr_funnel_plot(mr_results[[paste0(which_mr_run, "_mr_singlesnp")]])
        }
        
        if ("mr_leaveoneout" %in% mr_results[["run_info"]][["additional_analyses"]]) {
            plots[[paste0(which_mr_run, "_mr_leaveoneout")]] <- mr_leaveoneout_plot(mr_results[[paste0(which_mr_run, "_mr_leaveoneout")]])
        }
    }
    
    
    for (plot_type in names(plots)) {
        
        for (p in names(plots[[plot_type]])) {
            
            if (print_plots) {
                print(plots[[plot_type]][[p]])
            }
            
            if (!is.null(plot_dir)) {
                ggsave(paste0(plot_dir, "/", plot_type, "_", p, ".pdf"), plots[[plot_type]][[p]])
            }
        }
    }
    
    return(plots)
}

summarise_qtl_workflow <- function(
    run_settings,
    subfolder
) {
    # collect results 
    collated_results_df <- list()
    collated_results_df[[subfolder]] <- list()
    
    # which proteins were used in the run?
    protein_details <- read.csv(paste0("scripts/", run_settings[["rund"]], "/protein_files.tsv"), sep = "\t")
    
    # for each protein
    for (i in 1:nrow(protein_details)) {
        
        # this is where the protein's results are stored
        results_dir <- paste0("results/", run_settings[["rund"]], protein_details[i,1], "/")
        
        # if no pQTL results, move on...
        if (!file.exists(paste0(results_dir, subfolder, "/results.rds"))) {
            next
        }
        
        # get pQTL results
        qtl_res <- readRDS(paste0(results_dir, subfolder, "/results.rds"))

        for (which_mr_run in qtl_res[["run_info"]][["which_mr_runs"]]) {  # for cis, trans all...
            
            for (analysis_type in qtl_res[["run_info"]][["additional_analyses"]]) {  # for each mr method

                if (analysis_type == "ivw_radial") {
                    next
                }
                
                if(!is.null(qtl_res[[paste0(which_mr_run, "_", analysis_type)]])) {  # ensure run is not null
                    if (nrow(qtl_res[[paste0(which_mr_run, "_", analysis_type)]]) > 0) {  # ensure there are rows

                        # add result to results
                        to_add <- qtl_res[[paste0(which_mr_run, "_", analysis_type)]]
                        
                        if (!("package" %in% colnames(to_add))) {
                            to_add$package <- "TwoSampleMR"
                        }
                        
                        collated_results_df[[subfolder]][[analysis_type]] <- rbind(collated_results_df[[subfolder]][[analysis_type]], to_add)
                    }
                }
                
                if(!is.null(qtl_res[[paste0(which_mr_run, "_", analysis_type, "_cor")]])) {  # ensure run is not null
                    if (nrow(qtl_res[[paste0(which_mr_run, "_", analysis_type, "_cor")]]) > 0) {  # ensure there are rows

                        # add result to results
                        to_add <- qtl_res[[paste0(which_mr_run, "_", analysis_type)]]
                        
                        if (!("package" %in% colnames(to_add))) {
                            to_add$package <- "MendelianRandomization"
                        }
                        
                        collated_results_df[[subfolder]][[analysis_type]] <- rbind(collated_results_df[[subfolder]][[analysis_type]], to_add)
                    }
                }
            }

            if(!is.null(qtl_res[[paste0(which_mr_run, "_harm_removed_outliers")]])) {
                collated_results_df[[subfolder]][["harm_removed_outliers"]] <- rbind(collated_results_df[[subfolder]][["harm_removed_outliers"]],
                                                                                     qtl_res[[paste0(which_mr_run, "_harm_removed_outliers")]])
                
                collated_results_df[[subfolder]][[paste0("num_removed_outliers")]] <- rbind(collated_results_df[[subfolder]][[paste0("num_removed_outliers")]],
                                                                                            data.frame(exposure = paste0(protein_details[i,1], "_", which_mr_run),
                                                                                                       num_removed_outliers = qtl_res[["run_info"]][[paste0(which_mr_run, "_num_removed_outliers")]]))
            }
        }
    }
    
    # write the results (separate file for e and p QTLs)
    results_dir <- paste0("results/", run_settings[["rund"]], "/")
    
    for (analysis_type in c(qtl_res[["run_info"]][["additional_analyses"]], "harm_removed_outliers", "num_removed_outliers")) {
        
        if (analysis_type == "ivw_radial") {
            next
        }
        
        write.table(collated_results_df[[subfolder]][[analysis_type]], 
                    file = paste0("results/", run_settings[["rund"]], "/collated_", subfolder, "_", analysis_type, ".tsv"), 
                    row.names = F, quote = F, sep = "\t")
    }
    
    all_results <- list()
    for (package in c("TwoSampleMR", "MendelianRandomization")) {
        all_results[[package]] <- write_collated_df(results_dir, paste0("/collated_", subfolder, "_"), run_settings[["num_proteins"]], package)
    }
    
    return(all_results)
}

# collate results and write to file in results_dir
write_collated_df <- function(results_dir, results_file, num_proteins, package = "TwoSampleMR") {
    
    collated_tsvs <- paste0(results_dir, results_file)

    all_collated_mr <- read.csv(paste0(collated_tsvs, "mr.tsv"), sep = "\t")
    all_collated_mr$GeneID <- sapply(all_collated_mr$exposure, function(r) {strsplit(r, "_")[[1]][1]})
    all_collated_mr <- all_collated_mr[all_collated_mr$package == package | all_collated_mr$method == "Wald ratio",]
    
    collated_mr_heterogeneity <- read.csv(paste0(collated_tsvs, "mr_heterogeneity.tsv"), sep = "\t")
    collated_mr_heterogeneity$GeneID <- sapply(collated_mr_heterogeneity$exposure, function(r) {strsplit(r, "_")[[1]][1]})
    collated_mr_heterogeneity <- collated_mr_heterogeneity[collated_mr_heterogeneity$method == "Inverse variance weighted",]
    collated_mr_heterogeneity <- dplyr::select(collated_mr_heterogeneity, -method)
    
    collated_pleiotropy_test <- read.csv(paste0(collated_tsvs, "pleiotropy_test.tsv"), sep = "\t")
    colnames(collated_pleiotropy_test) <- c("id.exposure", "id.outcome", "outcome", "exposure", "egger_intercept", "egger_se", "egger_pval", "package")
    collated_pleiotropy_test$GeneID <- sapply(collated_pleiotropy_test$exposure, function(r) {strsplit(r, "_")[[1]][1]})
    
    collated_mr <- dplyr::left_join(all_collated_mr, collated_pleiotropy_test[collated_pleiotropy_test$package == package,])
    collated_mr <- dplyr::left_join(collated_mr, collated_mr_heterogeneity[collated_mr_heterogeneity$package == package,])

    if (package == "TwoSampleMR") {
        
        collated_directionality_test <- read.csv(paste0(collated_tsvs, "directionality_test.tsv"), sep = "\t")
        collated_directionality_test$GeneID <- sapply(collated_directionality_test$exposure, function(r) {strsplit(r, "_")[[1]][1]})
        
        collated_num_removed_outliers <- read.csv(paste0(collated_tsvs, "num_removed_outliers.tsv"), sep = "\t")
        collated_num_removed_outliers$GeneID <- sapply(collated_num_removed_outliers$exposure, function(r) {strsplit(r, "_")[[1]][1]})
        
        collated_mr <- dplyr::left_join(collated_mr, collated_directionality_test)
        collated_mr <- dplyr::left_join(collated_mr, collated_num_removed_outliers)
    }
    
    collated_mr <- collated_mr[collated_mr$method %in% c("Inverse variance weighted", "Inverse variance weighted (fixed effects)", "Wald ratio"),]
    collated_mr <- collated_mr[!(collated_mr$nsnp == 1 & collated_mr$method != "Wald ratio"),]
    collated_mr <- collated_mr[!(collated_mr$Q_pval > 0.05 & collated_mr$method == "Inverse variance weighted"),]
    collated_mr <- collated_mr[!(collated_mr$Q_pval <= 0.05 & collated_mr$method == "Inverse variance weighted (fixed effects)"),]
    
    collated_all <- collated_mr[grepl("all", collated_mr$exposure),]
    collated_cis <- collated_mr[grepl("cis", collated_mr$exposure),]
    
    collated_all$pval_fdr <- p.adjust(collated_all$pval, method = "fdr")
    collated_all$pval_bonf <- p.adjust(collated_all$pval, method = "bonferroni", n = num_proteins)
    collated_cis$pval_fdr <- p.adjust(collated_cis$pval, method = "fdr")
    
    if (nrow(collated_cis) > 0) {
        collated_cis$pval_bonf <- p.adjust(collated_cis$pval, method = "bonferroni", n = num_proteins)
    }
    
    for (method in c("Inverse variance weighted", "Inverse variance weighted (fixed effects)", "MR Egger", "Weighted median", "Weighted mode", "Maximum likelihood", "Robust adjusted profile score (RAPS)")) {
        collated_all <- add_method(collated_all, all_collated_mr, "all", method)
        collated_cis <- add_method(collated_cis, all_collated_mr, "cis", method)
    }

    stopifnot(length(unique(collated_all$GeneID)) == nrow(collated_all))
    collated_all <- collated_all[,which(apply(collated_all, 2, function(x) {!all(is.na(x))}))]
    write.csv(collated_all, paste0(results_dir, results_file, "_", package, "_all.csv"), row.names = FALSE)
    print(nrow(collated_all))
    
    stopifnot(length(unique(collated_cis$GeneID)) == nrow(collated_cis))
    collated_cis <- collated_cis[,which(apply(collated_cis, 2, function(x) {!all(is.na(x))}))]
    write.csv(collated_cis, paste0(results_dir, results_file, "_", package, "_cis.csv"), row.names = FALSE)
    print(nrow(collated_cis))
    
    return(list("all" = collated_all, "cis" = collated_cis))
}

# add method as a column
add_method <- function(mr_results, full_mr_results, qtl_type, method) {
    
    results_to_add <- full_mr_results[full_mr_results$method == method & grepl(qtl_type, full_mr_results$exposure),]
    
    results_to_add <- dplyr::select(results_to_add, "id.exposure", "id.outcome", "outcome", "GeneID", "nsnp", "b", "se", "pval")
    
    colnames(results_to_add) <- c("id.exposure", "id.outcome", "outcome", "GeneID", "nsnp", 
                                  paste0("b_", method), 
                                  paste0("se_", method), 
                                  paste0("pval_", method))
    
    mr_results <- dplyr::left_join(mr_results, results_to_add)
    
    return(mr_results)
}

ld_cor_mr <- function(mr_results, exposure_type, linkage_file, sink_plink, plink_bin) {
    
    if (nrow(mr_results[[paste0(exposure_type, "_harm")]]) <= 1) {
        
        return(mr_results)
        
    } else if (is.null(plink_bin)) {
        
        get_correlations <- correl <- FALSE
        
    } else {
        
        get_correlations <- correl <- TRUE
    }
    
    mr_results[[paste0(exposure_type, "_mr_object")]] <- dat_to_MRInput_modified(mr_results[[paste0(exposure_type, "_harm")]], 
                                                                                 get_correlations = get_correlations, 
                                                                                 linkage_file = linkage_file, 
                                                                                 sink_plink = sink_plink, 
                                                                                 plink_bin = plink_bin)[[1]]
    
    # mr ivw random
    mr_results[[paste0(exposure_type, "_ivw_cor")]] <- MendelianRandomization::mr_ivw(
        mr_results[[paste0(exposure_type, "_mr_object")]],
        model = "random",
        robust = FALSE,
        penalized = FALSE,
        correl = correl,
        weights = "simple"
    )
    
    mr_results[[paste0(exposure_type, "_mr")]] <- rbind(mr_results[[paste0(exposure_type, "_mr")]], data.frame(
        id.exposure = mr_results[[paste0(exposure_type, "_harm")]]$id.exposure[1], 
        id.outcome = mr_results[[paste0(exposure_type, "_harm")]]$id.outcome[1], 
        outcome = mr_results[[paste0(exposure_type, "_mr_object")]]@outcome, 
        exposure = mr_results[[paste0(exposure_type, "_mr_object")]]@exposure, 
        method = "Inverse variance weighted", 
        nsnp = length(mr_results[[paste0(exposure_type, "_mr_object")]]@snps), 
        b = mr_results[[paste0(exposure_type, "_ivw_cor")]]@Estimate, 
        se = mr_results[[paste0(exposure_type, "_ivw_cor")]]@StdError, 
        pval = mr_results[[paste0(exposure_type, "_ivw_cor")]]@Pvalue,
        package = "MendelianRandomization"
    ))
    
    # mr ivw fixed
    mr_results[[paste0(exposure_type, "_ivw_fe_cor")]] <- MendelianRandomization::mr_ivw(
        mr_results[[paste0(exposure_type, "_mr_object")]],
        model = "fixed",
        robust = FALSE,
        penalized = FALSE,
        correl = correl,
        weights = "simple"
    )
    
    mr_results[[paste0(exposure_type, "_mr")]] <- rbind(mr_results[[paste0(exposure_type, "_mr")]], data.frame(
        id.exposure = mr_results[[paste0(exposure_type, "_harm")]]$id.exposure[1], 
        id.outcome = mr_results[[paste0(exposure_type, "_harm")]]$id.outcome[1], 
        outcome = mr_results[[paste0(exposure_type, "_mr_object")]]@outcome, 
        exposure = mr_results[[paste0(exposure_type, "_mr_object")]]@exposure, 
        method = "Inverse variance weighted (fixed effects)", 
        nsnp = length(mr_results[[paste0(exposure_type, "_mr_object")]]@snps), 
        b = mr_results[[paste0(exposure_type, "_ivw_fe_cor")]]@Estimate, 
        se = mr_results[[paste0(exposure_type, "_ivw_fe_cor")]]@StdError, 
        pval = mr_results[[paste0(exposure_type, "_ivw_fe_cor")]]@Pvalue,
        package = "MendelianRandomization"
    ))
    
    # mr_maxlik
    mr_results[[paste0(exposure_type, "_mr_two_sample_ml_cor")]] <- MendelianRandomization::mr_maxlik(
        mr_results[[paste0(exposure_type, "_mr_object")]], 
        correl = correl
    )
    
    mr_results[[paste0(exposure_type, "_mr")]] <- rbind(mr_results[[paste0(exposure_type, "_mr")]], data.frame(
        id.exposure = mr_results[[paste0(exposure_type, "_harm")]]$id.exposure[1], 
        id.outcome = mr_results[[paste0(exposure_type, "_harm")]]$id.outcome[1], 
        outcome = mr_results[[paste0(exposure_type, "_mr_object")]]@outcome, 
        exposure = mr_results[[paste0(exposure_type, "_mr_object")]]@exposure, 
        method = "Maximum likelihood", 
        nsnp = length(mr_results[[paste0(exposure_type, "_mr_object")]]@snps), 
        b = mr_results[[paste0(exposure_type, "_mr_two_sample_ml_cor")]]@Estimate,
        se = mr_results[[paste0(exposure_type, "_mr_two_sample_ml_cor")]]@StdError, 
        pval = mr_results[[paste0(exposure_type, "_mr_two_sample_ml_cor")]]@Pvalue,
        package = "MendelianRandomization"
    ))
    
    if (length(mr_results[[paste0(exposure_type, "_mr_object")]]@snps) > 2) {
        
        # egger
        mr_results[[paste0(exposure_type, "_mr_egger_regression_cor")]] <- MendelianRandomization::mr_egger(
            mr_results[[paste0(exposure_type, "_mr_object")]], 
            robust = FALSE,
            penalized = FALSE,
            correl = correl
        )
        
        mr_results[[paste0(exposure_type, "_mr")]] <- rbind(mr_results[[paste0(exposure_type, "_mr")]], data.frame(
            id.exposure = mr_results[[paste0(exposure_type, "_harm")]]$id.exposure[1], 
            id.outcome = mr_results[[paste0(exposure_type, "_harm")]]$id.outcome[1], 
            outcome = mr_results[[paste0(exposure_type, "_mr_object")]]@outcome, 
            exposure = mr_results[[paste0(exposure_type, "_mr_object")]]@exposure, 
            method = "MR Egger", 
            nsnp = length(mr_results[[paste0(exposure_type, "_mr_object")]]@snps), 
            b = mr_results[[paste0(exposure_type, "_mr_egger_regression_cor")]]@Estimate,
            se = mr_results[[paste0(exposure_type, "_mr_egger_regression_cor")]]@StdError.Est, 
            pval = mr_results[[paste0(exposure_type, "_mr_egger_regression_cor")]]@Pvalue.Est,
            package = "MendelianRandomization"
        ))
        
        # mr_median
        mr_results[[paste0(exposure_type, "_mr_weighted_median_cor")]] <- MendelianRandomization::mr_median(
            mr_results[[paste0(exposure_type, "_mr_object")]], 
            weighting = "weighted"
        )
        
        mr_results[[paste0(exposure_type, "_mr")]] <- rbind(mr_results[[paste0(exposure_type, "_mr")]], data.frame(
            id.exposure = mr_results[[paste0(exposure_type, "_harm")]]$id.exposure[1], 
            id.outcome = mr_results[[paste0(exposure_type, "_harm")]]$id.outcome[1], 
            outcome = mr_results[[paste0(exposure_type, "_mr_object")]]@outcome, 
            exposure = mr_results[[paste0(exposure_type, "_mr_object")]]@exposure, 
            method = "Weighted median", 
            nsnp = length(mr_results[[paste0(exposure_type, "_mr_object")]]@snps), 
            b = mr_results[[paste0(exposure_type, "_mr_weighted_median_cor")]]@Estimate,
            se = mr_results[[paste0(exposure_type, "_mr_weighted_median_cor")]]@StdError, 
            pval = mr_results[[paste0(exposure_type, "_mr_weighted_median_cor")]]@Pvalue,
            package = "MendelianRandomization"
        ))
        
        # mr_mbe
        mr_results[[paste0(exposure_type, "_mr_weighted_mode_cor")]] <- MendelianRandomization::mr_median(
            mr_results[[paste0(exposure_type, "_mr_object")]], 
            weighting = "weighted"
        )
        
        mr_results[[paste0(exposure_type, "_mr")]] <- rbind(mr_results[[paste0(exposure_type, "_mr")]], data.frame(
            id.exposure = mr_results[[paste0(exposure_type, "_harm")]]$id.exposure[1], 
            id.outcome = mr_results[[paste0(exposure_type, "_harm")]]$id.outcome[1], 
            outcome = mr_results[[paste0(exposure_type, "_mr_object")]]@outcome, 
            exposure = mr_results[[paste0(exposure_type, "_mr_object")]]@exposure, 
            method = "Weighted mode", 
            nsnp = length(mr_results[[paste0(exposure_type, "_mr_object")]]@snps), 
            b = mr_results[[paste0(exposure_type, "_mr_weighted_mode_cor")]]@Estimate,
            se = mr_results[[paste0(exposure_type, "_mr_weighted_mode_cor")]]@StdError, 
            pval = mr_results[[paste0(exposure_type, "_mr_weighted_mode_cor")]]@Pvalue,
            package = "MendelianRandomization"
        ))
    }
    
    return(mr_results)
}
