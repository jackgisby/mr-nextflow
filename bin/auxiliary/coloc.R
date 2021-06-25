get_coloc_list <- function(gwas, list_type, LD = NULL, s = "", sdy="") {
    
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
    
    if (s != "") {
        coloc_in[["s"]] = as.numeric(s)
    }

    if (sdy != "") {
        coloc_in[["sdY"]] = as.numeric(sdy)
    }
    
    if (!is.null(LD)) {
        coloc_in[["LD"]] = LD
    }
    
    return(coloc_in)
}

process_sensitivity_data <- function(sensitivity_object, exposure_name, to, from, nsnps) {
    
    sensitivity_object$exposure <- exposure_name
    sensitivity_object$to <- to
    sensitivity_object$from <- from
    sensitivity_object$region_width <- to - from
    # sensitivity_object$p1 = 1e-4
    # sensitivity_object$p2 = 1e-4

    if (!("nsnps" %in% colnames(sensitivity_object))) {
        sensitivity_object$nsnps <- nsnps
    }
    
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
