apply_coloc <- function(
    exposure_name,
    exposure_gwas,
    outcome_gwas,
    results_dir,
    from,
    to,
    chr,
    linkage_files,
    plink_bin,
    make_plots = TRUE,
    p12=1e-5,
    s = 0.2,
    nref = 1000,
    susie = FALSE
) {
    colnames(outcome_gwas) <- gsub(".outcome", "", colnames(outcome_gwas))
    colnames(exposure_gwas) <- gsub(".exposure", "", colnames(exposure_gwas))
    
    # convert from TwoSampleMR format to coloc format
    exposure_coloc_in <- get_coloc_list(exposure_gwas, list_type = "quant")
    outcome_coloc_in <- get_coloc_list(outcome_gwas, list_type = "cc", s = s)
    
    # run coloc with default parameters
    coloc_res = coloc.abf(dataset1 = exposure_coloc_in, dataset2 = outcome_coloc_in, 
                          p1 = 1e-4, p2 = 1e-4, p12 = p12)
    
    # save rds
    saveRDS(coloc_res, paste0(results_dir, "/coloc_res.rds"))
    print(coloc_res)
    pdf(paste0(results_dir, "/sensitivity.pdf"))
    sensitivity_single <- sensitivity(coloc_res, rule = "H4 > 0.5")
    dev.off()
    
    sensitivity_single <- process_sensitivity_data(sensitivity_single, exposure_name, to, from, 1e-4, nrow(exposure_gwas))
    write.csv(sensitivity_single, paste0(results_dir, "/sensitivity.csv"), row.names = FALSE)
    
    if (!is.null(linkage_files) & susie) {
        
        # https://stephenslab.github.io/susieR/articles/finemapping_summary_statistics.html#fine-mapping-with-susier-using-summary-statistics-1
        
        LD = ld_matrix_modified(outcome_coloc_in$snp, chr, linkage_file=linkage_files, plink_bin=plink_bin, with_alleles = TRUE)
        
        exposure_coloc_in <- get_coloc_list(exposure_gwas, list_type = "quant", LD = LD)
        exposure_run_susie <- runsusie(exposure_coloc_in, nref = nref, p = 1e-4)
        saveRDS(exposure_run_susie, paste0(results_dir, "/exposure_susie.rds"))
        
        outcome_coloc_in <- get_coloc_list(outcome_gwas, list_type = "cc", LD = LD, s = s)
        outcome_run_susie <- runsusie(outcome_coloc_in, nref = nref, p = 1e-4)
        saveRDS(outcome_run_susie, paste0(results_dir, "/outcome_susie.rds"))
        
        susie_coloc_res <- coloc.susie(exposure_run_susie, outcome_run_susie, p1 = 1e-4, p2 = 1e-4, p12 = p12)
        saveRDS(multi_coloc_res, paste0(results_dir, "/coloc_susie.rds"))
        
        pdf(paste0(results_dir, "/sensitivity_susie.pdf"))
        susie_sensitivity <- sensitivity(susie_coloc_res, rule = "H4 > 0.5", dataset1 = exposure_coloc_in, dataset2 = outcome_coloc_in)
        dev.off()
        
        susie_sensitivity <- process_sensitivity_data(susie_sensitivity, exposure_name, to, from, 1e-4, nrow(exposure_gwas))
        write.csv(susie_sensitivity, paste0(results_dir, "/coloc_susie.csv"))
    }
    
    # ------------------------- plot
    
    if (make_plots & !is.null(linkage_files)) {
        
        try(make_variant_plot(
            exposure_name, 
            exposure_gwas,
            outcome_gwas,
            results_dir = results_dir,
            region_start = from,
            region_end = to,
            chr = chr,
            linkage_file = linkage_files[[chr]],
            plink_bin = plink_bin
        ))
    }
    
    return(sensitivity_single)
}

process_sensitivity_data <- function(sensitivity_object, exposure_name, to, from, p, nsnps) {
    
    sensitivity_object$exposure <- exposure_name
    sensitivity_object$region_width <- to - from
    sensitivity_object$p1 = 1e-4
    sensitivity_object$p2 = 1e-4
    sensitivity_object$nsnps <- nrow(exposure_gwas)
    
    # # p4 p3 ratio and prior for this
    # sensitivity_object$P4_P3_ratio <- sensitivity_object$PP.H4.abf / sensitivity_object$PP.H3.abf
    # sensitivity_object$P4_P3_prior_ratio <- ((10^8) * sensitivity_object$p12) / (sensitivity_object$nsnps - 1)
    
    # # priors for each hypothesis
    # sensitivity_object$p0 <- 1 - (sensitivity_object$p1 + sensitivity_object$p2 + sensitivity_object$p12)
    # sensitivity_object$p.H0 <-  sensitivity_object$p0 ^ sensitivity_object$nsnps
    # sensitivity_object$p.H1 <- sensitivity_object$nsnps * (sensitivity_object$p0 ^ (sensitivity_object$nsnps - 1)) * sensitivity_object$p1
    # sensitivity_object$p.H2 <- sensitivity_object$nsnps * (sensitivity_object$p0 ^ (sensitivity_object$nsnps - 1)) * sensitivity_object$p2
    # sensitivity_object$p.H3 <- sensitivity_object$nsnps * (sensitivity_object$nsnps - 1) * (sensitivity_object$p0 ^ (sensitivity_object$nsnps - 2)) * sensitivity_object$p1 * sensitivity_object$p2
    # sensitivity_object$p.H4 <- sensitivity_object$nsnps * (sensitivity_object$p0 ^ (sensitivity_object$nsnps - 1)) * sensitivity_object$p12

    return(sensitivity_object)
}

summarise_coloc_results <- function(
    run_settings,
    subfolder = "coloc"
) {
    # collect results 
    collated_coloc_res <- data.frame()
    collated_susie_res <- data.frame()
    
    # which proteins were used in the run?
    protein_details <- read.csv(paste0("scripts/", run_settings[["rund"]], "/protein_files.tsv"), sep = "\t")
    
    # for each protein
    for (i in 1:nrow(protein_details)) {
        
        # this is where the protein's results are stored
        results_dir <- paste0("results/", run_settings[["rund"]], protein_details[i,1], "/")

        # if no pQTL results, move on...
        if (!file.exists(paste0(results_dir, subfolder, "/sensitivity.csv"))) {
            next
        }
        
        # get pQTL results
        coloc_res <- read.csv(paste0(results_dir, subfolder, "/sensitivity.csv"))
        
        collated_coloc_res <- rbind(collated_coloc_res, coloc_res)
        
        # do same for susie
        if (!file.exists(paste0(results_dir, subfolder, "/coloc_susie.csv"))) {
            next
        }
        
        # get pQTL results
        susie_res <- read.csv(paste0(results_dir, subfolder, "/coloc_susie.csv"))
        
        collated_susie_res <- rbind(collated_susie_res, susie_res)
    }
    
    # write the results
    write.csv(collated_coloc_res, paste0("results/", run_settings[["rund"]], "/collated_", subfolder, ".csv"), row.names = FALSE)
    write.csv(collated_susie_res, paste0("results/", run_settings[["rund"]], "/collated_susie_", subfolder, ".csv"), row.names = FALSE)
    
    return(collated_coloc_res)
}

get_coloc_list <- function(gwas, list_type, LD = NULL, s = NULL) {
    
    if(!is.null(LD)) {
        stopifnot(colnames(LD) == rownames(LD))
        
        SNPs <- sapply(colnames(LD), function(SNP) {strsplit(SNP, "_")[[1]][1]})
        effect <- sapply(colnames(LD), function(SNP) {strsplit(SNP, "_")[[1]][2]})
        other <- sapply(colnames(LD), function(SNP) {strsplit(SNP, "_")[[1]][3]})
        
        if (!is.null(SNPs)) {
            gwas <- gwas[which(gwas$SNP %in% SNPs),]
        }
        
        to_remove <- vector("logical", length(SNPs))
        correct_orientation <- vector("logical", length(SNPs))
        for (i in 1:length(SNPs)) {
            
            snp_row <- gwas[gwas$SNP == SNPs[i],]
            stopifnot(nrow(snp_row) == 1)
            
            if (!(snp_row$effect_allele %in% c(effect[i], other[i]) | snp_row$other_allele %in% c(effect[i], other[i]))) {
                
                to_remove[i] <- TRUE
                
            } else if (snp_row$effect_allele == effect[i]) {
                
                correct_orientation[i] <- TRUE
                stopifnot(colnames(LD)[i] == paste(snp_row$SNP, snp_row$effect_allele, snp_row$other_allele, sep="_"))
                
            } else if (snp_row$effect_allele == other[i]) {
                
                correct_orientation[i] <- FALSE
                stopifnot(colnames(LD)[i] == paste(snp_row$SNP, snp_row$other_allele, snp_row$effect_allele, sep="_"))
                
            } else {
                stopifnot(FALSE)
            }
        }
        
        for (i in 1:length(SNPs)) {
            for (j in 1:length(SNPs)) {
                
                if (correct_orientation[i] != correct_orientation[j]) {
                    LD[i,j] <- - LD[i,j]
                }
            }
        }
        
        colnames(LD) <- rownames(LD) <- SNPs
        
        if (any(to_remove)) {
            message(paste("SNPs being removed due to incorrect alleles", SNPs[to_remove], sep = "\n"))
        }
        
        LD <- LD[which(!to_remove),which(!to_remove)]
    }
    
    coloc_in <- list(
        "pvalues" = gwas$pval,
        "N" = gwas$samplesize,
        "MAF" = gwas$eaf,
        "beta" = gwas$beta,
        "varbeta" = gwas$se ^ 2,
        "type" = list_type,
        "snp" = gwas$SNP,
        "position" = gwas$pos
    )
    
    if (!is.null(s)) {
        coloc_in[["s"]] = s
    }
    
    if (!is.null(LD)) {
        coloc_in[["LD"]] = LD
    }
    
    return(coloc_in)
}
