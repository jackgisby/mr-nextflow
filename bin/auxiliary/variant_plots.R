#' Generate position vs -log10(p) variant plots
#' 
#' Generates variant plots and locuszoom-like plots
#' 
#' @param exposure_name
#' Name of the exposure - e.g. gene name corresponding to a QTL.
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
#' If this parameter is NULL, does not plot locuszoom-like plots. Else, should be
#' file(s) for the calculation of linkage disequilibrium. If this is a single file,
#' this should be a string referring to the path of the file. If there is one file
#' per chromosome, then this should be a list; the keys of the list should be
#' integers referring to the chromosome number and the values should be a string
#' referring to the path of the file.
#' 
#' @param results_dir
#' Directory in which to write plot PDFs.
#' 
#' @param region_start
#' The minimum base at which to start plotting variants.
#' 
#' @param region_end
#' The maximum base at which to start plotting variants.
#' 
#' @param chr
#' The name of the chromosome to plot variants.
#' 
#' @param signif_line
#' A horizontal line will be drawn at `yintercept=-log10p(signif_line)`.
#' 
#' @return 
#' List of plots; keys correspond to the plot name and values correspond to
#' ggplot2 objects. 

make_variant_plot <- function(
    exposure_name,
    exposure_gwas,
    outcome_gwas,
    linkage_file = NULL,
    results_dir,
    region_start,
    region_end,
    chr,
    signif_line = 5e-8,
    print_plots = FALSE,
    plink_bin = NULL
) {
    exposure_gwas <- exposure_gwas[exposure_gwas$chr == chr &
                                   exposure_gwas$pos > region_start &
                                   exposure_gwas$pos < region_end,]
    
    if (nrow(exposure_gwas) == 0) {
        return()
    }
    
    exposure_gwas$mr_var <- "exposure"
    
    outcome_gwas <- outcome_gwas[outcome_gwas$chr == chr &
                                 outcome_gwas$pos > region_start &
                                 outcome_gwas$pos < region_end,]
    
    if (nrow(outcome_gwas) == 0) {
        return()
    }

    outcome_gwas$mr_var <- "outcome"

    plot_list <- list()
    
    if (!is.null(plink_bin)) {
        
        outcome_rsid <- outcome_gwas$SNP[which.min(outcome_gwas$pval)][1]
        exposure_rsid <- exposure_gwas$SNP[which.min(exposure_gwas$pval)][1]
        
        exposure_gwas$r2 <- get_lds(exposure_gwas, linkage_file, plink_bin, window_kb = region_end - region_start, rsid = exposure_rsid)
        outcome_gwas$r2 <- get_lds(outcome_gwas, linkage_file, plink_bin, window_kb = region_end - region_start, rsid = outcome_rsid)
        
        plot_list[["exposure_ld_plot"]] <- ggplot(exposure_gwas, aes(pos, -log10(pval), col = r2)) +
            geom_point(alpha = 0.7) +
            geom_hline(yintercept = -log10(signif_line), linetype = "dashed", alpha = 0.7) +
            ggtitle(paste0("ld: ", exposure_name, "; pval: ", exposure_name)) +
            scale_color_viridis() +
            geom_label_repel(data = exposure_gwas[exposure_gwas$SNP == exposure_rsid,], aes(label = exposure_rsid), col = "black", box.padding = 1, min.segment.length = unit(0, 'lines'))

        plot_list[["outcome_ld_plot"]] <- ggplot(outcome_gwas, aes(pos, -log10(pval), col = r2)) +
            geom_point(alpha = 0.7) +
            geom_hline(yintercept = -log10(signif_line), linetype = "dashed", alpha = 0.7) +
            ggtitle("ld: outcome; pval: outcome") +
            scale_color_viridis() +
            geom_label_repel(data = outcome_gwas[outcome_gwas$SNP == outcome_rsid,], aes(label = outcome_rsid), col = "black", box.padding = 1, min.segment.length = unit(0, 'lines'))
    }
    
    if (!is.null(plink_bin)) {

        exposure_gwas$outcome_r2 <- get_lds(exposure_gwas, linkage_file, plink_bin, 
                                            window_kb = region_end - region_start, 
                                            rsid = outcome_rsid)
        
        outcome_gwas$exposure_r2 <- get_lds(outcome_gwas, linkage_file, plink_bin, 
                                            window_kb = region_end - region_start, 
                                            rsid = exposure_rsid)
        
        plot_list[["exposure_overlayed_ld_plot"]] <- ggplot(exposure_gwas, aes(pos, -log10(pval), col = outcome_r2)) +
            geom_point(alpha = 0.7) +
            geom_hline(yintercept = -log10(signif_line), linetype = "dashed", alpha = 0.7) +
            ggtitle(paste0("ld: outcome; pval: ", exposure_name)) +
            scale_color_viridis() +
            geom_label_repel(data = exposure_gwas[exposure_gwas$SNP == outcome_rsid,], aes(label = outcome_rsid), col = "black", box.padding = 1, min.segment.length = unit(0, 'lines'))
        
        plot_list[["outcome_overlayed_ld_plot"]] <- ggplot(outcome_gwas, aes(pos, -log10(pval), col = exposure_r2)) +
            geom_point(alpha = 0.7) +
            geom_hline(yintercept = -log10(signif_line), linetype = "dashed", alpha = 0.7) +
            ggtitle(paste0("ld: ", exposure_name, "; pval: outcome")) +
            scale_color_viridis() +
            geom_label_repel(data = outcome_gwas[outcome_gwas$SNP == exposure_rsid,], aes(label = exposure_rsid), col = "black", box.padding = 1, min.segment.length = unit(0, 'lines'))
    }
    
    all_gwas <- dplyr::bind_rows(outcome_gwas, exposure_gwas)

    # make standard plot
    variant_plot <- ggplot(all_gwas, aes(pos, -log10(pval), col = mr_var)) +
        geom_point(alpha = 0.6) +
        geom_hline(yintercept = -log10(signif_line), linetype = "dashed", alpha = 0.7) +
        ggtitle(paste0("exposure: ", exposure_name))
    
    # make relative plot
    exposure_gwas$pval <- -log10(exposure_gwas$pval) / max(-log10(exposure_gwas$pval))
    outcome_gwas$pval <- -log10(outcome_gwas$pval) / max(-log10(outcome_gwas$pval))
    all_gwas <- dplyr::bind_rows(outcome_gwas, exposure_gwas)
    
    rel_variant_plot <- ggplot(all_gwas, aes(pos, pval, col = mr_var)) +
        geom_point(alpha = 0.6) +
        ggtitle(paste0("exposure: ", exposure_name)) +
        ylab("-log10(pval) / max(-log10(pval))")
    
    plot_list[["rel_variant_plot"]] <- rel_variant_plot
    plot_list[["variant_plot"]] <- variant_plot
    
    # save to results
    for (plot_name in names(plot_list)) {
        
        if (print_plots) {
            print(plot_list[[plot_name]])
        }
        
        ggsave(paste0(results_dir, "/", exposure_name, "_", plot_name, ".pdf"), plot_list[[plot_name]])
    }

    return(plot_list)
}

#' get lds for plotting
#' 
#' gets the ld with the reference variant (whichever has the minimal p value)
#' using plink and a supplied file
#' 
#' @param gwas
#' A `data.frame` in the format returned by \code{\link[TwoSampleMR]{format_data}}.
#' Contains the variants for the exposure, for instance a QTL, or outcome.
#' 
#' @param linkage_file
#' A BED file for use by plink to extract LD.
#' 
#' @param plink_bin
#' Location of the plink executable.
#' 
#' @param window_kb
#' The window around which to search for the variant (refers to a single side,
#' not the total window)
#' 
#' @return 
#' Returns the `gwas` `data.frame` with an additional column, r2, containing
#' the linkage correlation with the most significant variant
#' 

get_lds <- function(
    gwas, 
    linkage_file, 
    plink_bin,
    window_kb,
    rsid=NULL
) {
    
    if (is.null(rsid)) {
        rsid <- gwas$SNP[which.min(gwas$pval)][1]
    }
    
    r2 <- get_variant_r2(rsid, linkage_file, window_kb * 2, plink_bin)
    
    gwas <- dplyr::left_join(gwas, r2, by=c("SNP"="SNP_B"))
    
    return(gwas$R2)
}
