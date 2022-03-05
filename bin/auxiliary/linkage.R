#' Wrapper around plink-clump for TwoSampleMR data
#' 
#' @param to_prune
#' A `data.frame` in the format returned by \code{\link[TwoSampleMR]{format_data}}.
#' Contains the "significant" variants for the exposure, for instance a QTL.
#' 
#' @param plink_bin
#' Path to PLINK program. If NULL, uses an API to access data. Avoid repeated
#' calls - mainly for off-cluster debugging.
#' 
#' @param linkage_file
#' File(s) for the calculation of linkage disequilibrium. If this is a single file,
#' this should be a string referring to the path of the file. If there is one file
#' per chromosome, then this should be a list; the keys of the list should be
#' integers referring to the chromosome number and the values should be a string
#' referring to the path of the file. Set to NULL if plink_bin is NULL.
#' 
#' @param clump_kb
#' Parameter clump-kb for PLINK clumping. 
#' 
#' @param clump_r2
#' Parameter clump-r2 for PLINK clumping.
#' 
#' @param clump_p1
#' Parameter clump-p1 for PLINK clumping.
#' 
#' @param clump_p2
#' Parameter clump-p2 for PLINK clumping.
#' 
#' @param file_format
#' The file format to be used as input to PLINK.
#' 
#' @references 
#' Adapted from \code{\link[TwoSampleMR]{clump_data}}.
#' 
#' @return 
#' Returns the original `to_prune` `data.frame` with pruned variants.

plink_clump <- function(
    to_prune, 
    linkage_file,  # list of files for each chromosome or single file
    clump_kb = 10000, 
    clump_r2 = 0.01, 
    clump_p = 1, 
    plink_bin = "plink",
    plink_memory = 15000
) {
    if (plink_bin == "") {
        message("plink_bin is NULL: using API (pop = EUR)")
        linkage_file <- NULL
    }
    
    if (nrow(to_prune) == 0) {
        return()
    }
    
    pval_column <- "pval.exposure"
    
    if (!is.data.frame(to_prune)) {
        stop("Expecting data frame returned from format_data")
    }
    
    if ("pval.exposure" %in% names(to_prune) & "pval.outcome" %in% names(to_prune)) {
        
        message("pval.exposure and pval.outcome columns present. Using pval.exposure for clumping.")
        
    } else if (!"pval.exposure" %in% names(to_prune) & "pval.outcome" %in% names(to_prune)) {
        
        message("pval.exposure column not present, using pval.outcome column for clumping.")
        pval_column <- "pval.outcome"
        
    } else if (!"pval.exposure" %in% names(to_prune)) {
        
        message("pval.exposure not present, setting clumping p-value to 0.99 for all variants")
        to_prune$pval.exposure <- 0.99
        
    } else {
        pval_column <- "pval.exposure"
    }
    
    if (!"id.exposure" %in% names(to_prune)) {
        to_prune$id.exposure <- random_string(1)
    }
    
    d <- data.frame(rsid = to_prune$SNP, pval = to_prune[[pval_column]])

    if (typeof(linkage_file) == "list") {
        out = c()
        
        for (i in as.integer(unique(to_prune$chr))) {  # for each chr: use chr-specific file
            
            chr_out <- try(ld_clump_modified(
                dat = d,
                clump_kb = clump_kb, 
                clump_r2 = clump_r2, 
                clump_p = clump_p,
                bfile = linkage_file[[i]],
                plink_bin = plink_bin,
                plink_memory = plink_memory
            ))
            
            if (all(class(chr_out) != "try-error")) {
                out <- rbind(out, chr_out)
            }
        }

    } else if (typeof(linkage_file) == "character") {
        out <- try(ld_clump_modified(
            dat = d,
            clump_kb = clump_kb, 
            clump_r2 = clump_r2, 
            clump_p = clump_p,
            bfile = linkage_file,
            plink_bin = plink_bin,
            plink_memory = plink_memory
        ))

    } else {
        out <- try(ieugwasr::ld_clump(
            dat = d,
            clump_kb = clump_kb, 
            clump_r2 = clump_r2, 
            clump_p = clump_p,
            bfile = linkage_file,
            plink_bin = plink_bin
        ))
    }
    
    keep <- to_prune$SNP %in% out$rsid
    
    return(to_prune[keep,])
}

ld_clump_modified <- function (dat = NULL, clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.99, 
                               pop = "EUR", access_token = NULL, bfile = NULL, plink_bin = NULL,
                               plink_memory = 15000) 
{
    stopifnot("rsid" %in% names(dat))
    stopifnot(is.data.frame(dat))

    if (is.null(bfile)) {
        message("Please look at vignettes for options on running this locally if you need to run many instances of this command.")
    }

    if (!"pval" %in% names(dat)) {
        if ("p" %in% names(dat)) {
            warning("No 'pval' column found in dat object. Using 'p' column.")
            dat[["pval"]] <- dat[["p"]]
        }
        else {
            warning("No 'pval' column found in dat object. Setting p-values for all SNPs to clump_p parameter.")
            dat[["pval"]] <- clump_p
        }
    }

    if (!"id" %in% names(dat)) {
        dat$id <- random_string(1)
    }

    if (is.null(bfile)) {
        access_token = check_access_token()
    }

    ids <- unique(dat[["id"]])
    res <- list()
    for (i in 1:length(ids)) {
        x <- subset(dat, dat[["id"]] == ids[i])
        if (nrow(x) == 1) {
            message("Only one SNP for ", ids[i])
            res[[i]] <- x
        }
        else {
            message("Clumping ", ids[i], ", ", nrow(x), " variants, using ", 
                    pop, " population reference")
            
            res[[i]] <- modified_ld_clump_local(x, clump_kb = clump_kb, 
                                                clump_r2 = clump_r2, clump_p = clump_p, bfile = bfile, 
                                                plink_bin = plink_bin, plink_memory = plink_memory)
        }
    }
    res <- dplyr::bind_rows(res)
    return(res)
}

modified_ld_clump_local <- function(dat, clump_kb, clump_r2, clump_p, bfile, plink_bin, plink_memory)
{
    
    # Make textfile
    shell <- ifelse(Sys.info()['sysname'] == "Windows", "cmd", "sh")
    fn <- tempfile()
    write.table(data.frame(SNP=dat[["rsid"]], P=dat[["pval"]]), file=fn, row.names=F, col.names=T, quote=F)
    
    fun2 <- paste0(
        shQuote(plink_bin, type=shell),
        " --bfile ", shQuote(bfile, type=shell),
        " --clump ", shQuote(fn, type=shell), 
        " --clump-p1 ", clump_p, 
        " --clump-r2 ", clump_r2, 
        " --clump-kb ", clump_kb, 
        " --memory ", plink_memory,  # add memory flag so that it doesn't eat up memory
        " --out ", shQuote(fn, type=shell)
    )
    system(fun2)

    res <- read.table(paste(fn, ".clumped", sep=""), header=T)
    unlink(paste(fn, "*", sep=""))
    y <- subset(dat, !dat[["rsid"]] %in% res[["SNP"]])

    if(nrow(y) > 0) {
        message("Removing ", length(y[["rsid"]]), " of ", nrow(dat), " variants due to LD with other variants or absence from LD reference panel")
    }

    return(subset(dat, dat[["rsid"]] %in% res[["SNP"]]))
}

random_string <- function(n=1, len=6)
{
    randomString <- c(1:n)
    for (i in 1:n)
    {
        randomString[i] <- paste(sample(c(0:9, letters, LETTERS),
                                        len, replace=TRUE),
                                 collapse="")
    }
    return(randomString)
}

#' Uses plink to get r2 to a reference variant
#' 
#' @param variant
#' rsid of the variant.
#' 
#' @param bfile
#' The BED file for plink to extract r2 statistics from.
#' 
#' @param window_kb
#' Plink's --ld-window-kb parameter
#' 
#' @param plink_bin
#' The location of the plink executable
#' 
#' @return 
#' A 'data.frame' containing rsid and r2 values generated by plink.

get_variant_r2 <- function(
    variant, 
    bfile, 
    window_kb, 
    plink_bin,
    plink_memory = 5000
) {
    shell <- ifelse(Sys.info()["sysname"] == "Windows", "cmd", "sh")
    fn <- tempfile()
    
    system(paste0(shQuote(plink_bin, type = shell),
                  " --bfile ", shQuote(bfile, type = shell),
                  " --r2",
                  " --ld-snp ", variant,
                  " --ld-window-kb ", window_kb,
                  " --ld-window 99999",
                  " --ld-window-r2 0",
                  " --memory ", plink_memory,
                  " --out ", shQuote(fn, type = shell)
    ))
    
    res <- read.table(paste(fn, ".ld", sep = ""), header = T)
    unlink(paste(fn, "*", sep = ""))
    
    return(res)
}

#' Get LD matrix for list of variants
#'
#' This function takes a list of variants and searches for them in samples from 1000 Genomes phase 3 data
#' It then creates an LD matrix of r values (signed, and not squared)
#' All LD values are with respect to the major alleles in the 1000G dataset. You can specify whether the allele names are displayed
#'
#' @param variants List of variants (rsids)
#' @param with_alleles Whether to append the allele names to the SNP names. Default: TRUE
#' @param pop Super-population to use as reference panel. Default = "EUR". Options are EUR, SAS, EAS, AFR, AMR. 'legacy' also available - which is a previously used verison of the EUR panel with a slightly different set of markers
#' @param bfile If this is provided then will use the API. Default = NULL
#' @param plink_bin If null and bfile is not null then will detect packaged plink binary for specific OS. Otherwise specify path to plink binary. Default = NULL
#'
#' @export
#' @return Matrix of LD r values

ld_matrix_modified <- function(
    variants, 
    chr,
    linkage_file=NULL,  # list of bfile locations for each chromosome
    with_alleles=TRUE, 
    pop="EUR", 
    plink_bin=NULL,
    plink_memory=5000
)  {
    if (length(variants) > 500 & is.null(linkage_file)) {
        stop("SNP list must be smaller than 500. Try running locally by providing local ld reference with bfile argument. See vignettes for a guide on how to do this.")
    }
    
    if (is.null(linkage_file)) {
        message("Please look at vignettes for options on running this locally if you need to run many instances of this command.")
    }
    
    if (!is.null(linkage_file)) {
        
        if (is.null(chr)) {  # chr given if linkage_file is a list of file locations - can use chr to select relevant file
            LD <- ld_matrix_local_modified(variants, bfile=linkage_file, plink_bin=plink_bin, with_alleles=with_alleles, plink_memory=plink_memory)
        } else {
            LD <- ld_matrix_local_modified(variants, bfile=linkage_file[[chr]], plink_bin=plink_bin, with_alleles=with_alleles, plink_memory=plink_memory)
        }
        
        return(LD)
    }
    
    res <- ieugwasr::get_query_content(ieugwasr::api_query('ld/matrix', query = list(rsid = variants, pop = pop), access_token=NULL))
    
    if (all(is.na(res))) {
        stop("None of the requested variants were found")
    }
    
    variants2 <- res$snplist
    res <- res$matrix
    res <- matrix(as.numeric(res), nrow(res), ncol(res))
    variants3 <- do.call(rbind, strsplit(variants2, split="_"))
    
    if (with_alleles) {
        rownames(res) <- variants2
        colnames(res) <- variants2
    } else {
        rownames(res) <- variants3[,1]
        colnames(res) <- variants3[,1]
    }
    
    missing <- variants[!variants %in% variants3[,1]]
    
    if (length(missing) > 0) {
        warning("The following variants are not present in the LD reference panel\n", paste(missing, collapse="\n"))
    }
    
    ord <- match(variants3[,1], variants)
    res <- res[order(ord), order(ord)]
    
    return(res)
}

multi_chr_ld_matrix <- function(
    x, 
    linkage_file=NULL,  # list of bfile locations for each chromosome
    plink_bin=NULL,
    plink_memory=5000
) {
    if (is.null(linkage_file)) {
        return(ld_matrix_modified(x$SNP, NULL, with_alleles = TRUE, linkage_file = linkage_file, plink_bin = plink_bin))
    }
    
    if (typeof(linkage_file) == "list") {  # list of linkage files, one for each chromosome

        LD <- matrix(nrow=0, ncol=0)

        for (chr in as.integer(unique(x$chr.exposure))) {
            LD_chr <- ld_matrix_modified(x$SNP[x$chr.exposure == chr], chr,
                                        linkage_file=linkage_file,  # list of bfile locations for each chromosome
                                        with_alleles=TRUE, plink_bin=plink_bin)
            
            if (ncol(LD) > 0) {
                
                new_names <- c(colnames(LD), colnames(LD_chr))
                
                LD_chr <- cbind(matrix(0, nrow=nrow(LD_chr), ncol=nrow(LD)), LD_chr)
                LD <- cbind(LD, matrix(0, nrow=nrow(LD), ncol=nrow(LD_chr)))
                LD <- rbind(LD, LD_chr)
                
                colnames(LD) <- new_names
                rownames(LD) <- new_names
                
            } else {
                LD <- LD_chr
            }
        }
    } else if (typeof(linkage_file) == "character") {
        LD <- ld_matrix_modified(x$SNP, NULL, linkage_file=linkage_file,  # list of bfile locations for each chromosome
                                 with_alleles=TRUE, plink_bin=plink_bin)
    }

    stopifnot(colnames(LD) == rownames(LD))

    return(LD)
}

#' Get LD matrix using local plink binary and reference dataset
#'
#' @param variants List of variants (rsids)
#' @param bfile Path to bed/bim/fam ld reference panel
#' @param plink_bin Specify path to plink binary. Default = NULL. See https://github.com/explodecomputer/plinkbinr for convenient access to plink binaries
#' @param with_alleles Whether to append the allele names to the SNP names. Default: TRUE
#'
#' @export
#' @return data frame
ld_matrix_local_modified <- function(variants, bfile, plink_bin, with_alleles=FALSE, plink_memory=5000) {
    
    # Make textfile
    shell <- ifelse(Sys.info()['sysname'] == "Windows", "cmd", "sh")
    fn <- tempfile()
    write.table(data.frame(variants), file=fn, row.names=F, col.names=F, quote=F)
    
    fun1 <- paste0(
        shQuote(plink_bin, type=shell),
        " --bfile ", shQuote(bfile, type=shell),
        " --extract ", shQuote(fn, type=shell), 
        " --make-just-bim ", 
        " --memory ", plink_memory,  # add memory flag so that it doesn't eat up memory, 
        " --out ", shQuote(fn, type=shell)
    )
    
    system(fun1)
    bim <- read.table(paste0(fn, ".bim"), stringsAsFactors=FALSE, colClasses = rep("character", 6))
    
    fun2 <- paste0(
        shQuote(plink_bin, type=shell),
        " --bfile ", shQuote(bfile, type=shell),
        " --extract ", shQuote(fn, type=shell), 
        " --r square ", 
        " --memory ", plink_memory,  # add memory flag so that it doesn't eat up memory, 
        " --out ", shQuote(fn, type=shell)
    )
    
    system(fun2)
    res <- as.matrix(read.table(paste0(fn, ".ld"), header=FALSE))
    
    if (with_alleles) {
        rownames(res) <- colnames(res) <- paste(bim$V2, bim$V5, bim$V6, sep="_")
    } else {
        rownames(res) <- colnames(res) <- bim$V2
    }
    
    return(res)
}

#' Convert TwoSampleMR format to MendelianRandomization format
#'
#' The MendelianRandomization package offers MR methods that can 
#' be used with the same data used in the TwoSampleMR package. This
#' function converts from the TwoSampleMR format to the MRInput class.
#'
#' @param dat Output from the [`harmonise_data`] function.
#' @param get_correlations Default `FALSE`. If `TRUE` then extract the LD matrix for the SNPs from the European 1000 genomes data on the MR-Base server.
#'
#' @export
#' @return List of MRInput objects for each exposure/outcome combination

dat_to_MRInput_modified <- function(
    dat, 
    get_correlations = TRUE,
    linkage_file, 
    plink_bin
)  {
    
    if (is.null(plink_bin)) {
        message("Plink bin is NULL, not getting LD correlations")
        
        linkage_file <- NULL
        get_correlations <- FALSE

    }
    
    out <- plyr::dlply(dat, c("exposure", "outcome"), function(x) {
        
        x <- plyr::mutate(x)
        message("Converting:")
        message(" - exposure: ", x$exposure[1])
        message(" - outcome: ", x$outcome[1])
        
        if (get_correlations) {
            
            message(" - obtaining LD matrix")
            
            ld <- multi_chr_ld_matrix(x, linkage_file=linkage_file, plink_bin=plink_bin)
            
            out <- harmonise_ld_dat(x, ld)
            
            if (is.null(out)) {
                return(NULL)
            }
            
            x <- out$x
            ld <- out$ld
            MendelianRandomization::mr_input(bx = x$beta.exposure, 
                                             bxse = x$se.exposure, by = x$beta.outcome, byse = x$se.outcome, 
                                             exposure = x$exposure[1], outcome = x$outcome[1], 
                                             snps = x$SNP, effect_allele = x$effect_allele.exposure, 
                                             other_allele = x$other_allele.exposure, eaf = x$eaf.exposure, 
                                             correlation = ld)
            
        } else {
            
            MendelianRandomization::mr_input(bx = x$beta.exposure, 
                                             bxse = x$se.exposure, by = x$beta.outcome, byse = x$se.outcome, 
                                             exposure = x$exposure[1], outcome = x$outcome[1], 
                                             snps = x$SNP, effect_allele = x$effect_allele.exposure, 
                                             other_allele = x$other_allele.exposure, eaf = x$eaf.exposure)
        }
    })
    
    return(out)
}

harmonise_ld_dat <- function(
    x,
    ld
)  {
    if (nrow(x) > 1) {

        snpnames <- do.call(rbind, strsplit(rownames(ld), split = "_"))

        # remove those in ld that are not also in x
        i1 <- snpnames[, 1] %in% x$SNP
        ld <- ld[i1, i1]
        
        snpnames <- do.call(rbind, strsplit(rownames(ld), split = "_"))
        
        # remove those in x that are not also in ld
        i2 <- x$SNP %in% snpnames[, 1]
        x <- x[i2, ]
        
        # reorder ld by x
        i3 <- match(x$SNP, snpnames[, 1])
        ld <- ld[i3, i3]
    }
    
    snpnames <- do.call(rbind, strsplit(rownames(ld), split = "_"))

    stopifnot(colnames(ld) == rownames(ld))
    stopifnot(all(snpnames[, 1] == x$SNP))
    x$effect_allele.exposure <- as.character(x$effect_allele.exposure)
    x$other_allele.exposure <- as.character(x$other_allele.exposure)
    
    # Set1 x and ld alleles match
    snpnames <- data.frame(snpnames, stringsAsFactors = FALSE)
    snpnames <- merge(subset(x, select = c(SNP, effect_allele.exposure, other_allele.exposure)), snpnames, by.x="SNP", by.y="X1")
    snpnames <- snpnames[match(x$SNP, snpnames$SNP) ,]
    
    snpnames$keep <- (snpnames$X2 == snpnames$effect_allele.exposure & snpnames$X3 == snpnames$other_allele.exposure) |
        (snpnames$X3 == snpnames$effect_allele.exposure & snpnames$X2 == snpnames$other_allele.exposure)
    
    # What happens if everything is gone?
    if (nrow(x) == 0) {
        message(" - none of the SNPs could be aligned to the LD reference panel")
        return(NULL)
    }
    
    if (any(!snpnames$keep)) {
        message(" - the following SNPs could not be aligned to the LD reference panel: \n- ", paste(subset(snpnames, !keep)$SNP, collapse="\n - "))
    }
    
    snpnames$flip1 <- snpnames$X2 != snpnames$effect_allele.exposure
    
    x <- subset(x, SNP %in% snpnames$SNP)
    
    temp1 <- x$effect_allele.exposure[snpnames$flip1]
    temp2 <- x$other_allele.exposure[snpnames$flip1]
    
    x$beta.exposure[snpnames$flip1] <- x$beta.exposure[snpnames$flip1] * -1
    x$beta.outcome[snpnames$flip1] <- x$beta.outcome[snpnames$flip1] * -1
    
    x$effect_allele.exposure[snpnames$flip1] <- temp2
    x$other_allele.exposure[snpnames$flip1] <- temp1
    
    rownames(ld) <- snpnames$SNP
    colnames(ld) <- snpnames$SNP
    
    if (any(!snpnames$keep)) {
        message("Removing ", sum(!snpnames$keep), " variants due to harmonisation issues")
        ld <- ld[snpnames$keep, snpnames$keep]
        x <- x[snpnames$keep, ]
    }
    
    return(list(x = x, ld = ld))
}
