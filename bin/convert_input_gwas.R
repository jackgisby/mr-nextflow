# load libraries
library("optparse")
library("data.table")
library("TwoSampleMR")

# do argument parsing
option_list = list(
    make_option(c("--exposure_input_data"), type="character", default=NULL, help="exposure input file name", metavar="character"),
    make_option(c("--outcome_input_data"), type="character", default=NULL, help="outcome input file name", metavar="character"),
    make_option(c("--auxiliary_script_dir"), type="character", default=NULL, help="the location of helper scripts", metavar="character"),
    make_option(c("--p_cutoff"), type="numeric", default=NULL),
    make_option(c("--cis_region"), type="numeric", default=NULL),
    make_option(c("--plink_memory"), type="numeric", default=NULL),
    make_option(c("--plink_clump_r2"), type="numeric", default=NULL),
    make_option(c("--plink_clump_kb"), type="numeric", default=NULL),
    make_option(c("--plink_bin"), type="character", default=NULL),
    make_option(c("--plink_linkage_files"), type="character", default=NULL),
    make_option(c("--gene_positions"), type="character", default=NULL),
    make_option(c("--gene_filenames"), type="character", default=NULL),
    make_option(c("--exposure_snp_col"), type="character", default=NULL),
    make_option(c("--exposure_beta_col"), type="character", default=NULL),
    make_option(c("--exposure_se_col"), type="character", default=NULL),
    make_option(c("--exposure_eaf_col"), type="character", default=NULL),
    make_option(c("--exposure_effect_allele_col"), type="character", default=NULL),
    make_option(c("--exposure_other_allele_col"), type="character", default=NULL),
    make_option(c("--exposure_pval_col"), type="character", default=NULL),
    make_option(c("--exposure_samplesize_col"), type="character", default=NULL),
    make_option(c("--exposure_chr_col"), type="character", default=NULL),
    make_option(c("--exposure_pos_col"), type="character", default=NULL),
    make_option(c("--outcome_snp_col"), type="character", default=NULL),
    make_option(c("--outcome_beta_col"), type="character", default=NULL),
    make_option(c("--outcome_se_col"), type="character", default=NULL),
    make_option(c("--outcome_eaf_col"), type="character", default=NULL),
    make_option(c("--outcome_effect_allele_col"), type="character", default=NULL),
    make_option(c("--outcome_other_allele_col"), type="character", default=NULL),
    make_option(c("--outcome_pval_col"), type="character", default=NULL),
    make_option(c("--outcome_samplesize_col"), type="character", default=NULL),
    make_option(c("--outcome_chr_col"), type="character", default=NULL),
    make_option(c("--outcome_pos_col"), type="character", default=NULL),
    make_option(c("--outcome_samplesize"), type="numeric", default=NULL)
)
 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

print("options")
print(opt)

source(paste0(opt$auxiliary_script_dir, "/linkage.R"))

convert_input_gwas <- function(opt) {

    # get gene name
    gene_to_filename <- data.frame(fread(opt$gene_filenames))
    gene_to_filename <- gene_to_filename[gene_to_filename$filename == opt$exposure_input_data,]

    stopifnot(nrow(gene_to_filename) == 1)
    exposure_name <- tolower(gene_to_filename$gene_id[1])

    cis_name <- paste(exposure_name, "cis", "harm", sep="_")
    all_name <- paste(exposure_name, "all", "harm", sep="_")
    blank_return <- list(cis_name=data.frame(), all_name=data.frame())
    names(blank_return) <- c(cis_name, all_name)

    # load exposure w/ fread
    exposure_gwas <- fread(opt$exposure_input_data)

    print("exposure")
    print(head(exposure_gwas))

    # make list to map exposure colnames
    exposure_cols = list("SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "pval", "samplesize", "chr_col", "pos")
    old_exposure_cols <- c(
        opt$exposure_snp_col, 
        opt$exposure_beta_col, 
        opt$exposure_se_col, 
        opt$exposure_eaf_col, 
        opt$exposure_effect_allele_col, 
        opt$exposure_other_allele_col, 
        opt$exposure_pval_col,
        opt$exposure_samplesize_col, 
        opt$exposure_chr_col, 
        opt$exposure_pos_col
    )
    names(exposure_cols) <- old_exposure_cols

    # select exposure columns
    exposure_gwas <- exposure_gwas[, ..old_exposure_cols]

    print("select exposure colnames")
    print(head(exposure_gwas))

    # modify exposure colnames
    for (old_name in names(exposure_cols)) {
        setnames(exposure_gwas, old_name, exposure_cols[[old_name]])
    }

    print("change exposure colnames")
    print(head(exposure_gwas))

    # remove invalid entries
    exposure_gwas <- exposure_gwas[exposure_gwas$SNP != ".",]
    exposure_gwas <- exposure_gwas[exposure_gwas$SNP != "",]
    exposure_gwas <- exposure_gwas[!is.na(exposure_gwas$pval),]

    print("remove invalid exposures")
    print(head(exposure_gwas))

    # limit exposure to significant variants
    exposure_gwas <- data.frame(exposure_gwas[exposure_gwas$pval <= opt$p_cutoff,])

    print("limited to significant exposures")
    print(head(exposure_gwas))

    if (length(exposure_gwas) == 0) {
        message("no significant exposure variants")
        return(list(cis_name=data.frame(), all_name=data.frame()))
    }

    # remove snps in MHC region
    which_mhc = exposure_gwas$chr_col == 6 & exposure_gwas$pos > 26000000 & exposure_gwas$pos < 34000000
    exposure_gwas <- exposure_gwas[-which_mhc,]

    if (length(exposure_gwas) == 0) {
        message("no significant exposure variants following removal of MHC")
        return()
    }

    # load outcome w/ fread
    outcome_gwas <- fread(opt$outcome_input_data)

    # make colname list using input options
    outcome_cols = list("SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "pval", "samplesize", "chr_col", "pos")
    old_outcome_cols <- c(
        opt$outcome_snp_col, 
        opt$outcome_beta_col, 
        opt$outcome_se_col, 
        opt$outcome_eaf_col, 
        opt$outcome_effect_allele_col, 
        opt$outcome_other_allele_col, 
        opt$outcome_pval_col,
        opt$outcome_samplesize_col, 
        opt$outcome_chr_col, 
        opt$outcome_pos_col
    )
    names(outcome_cols) <- old_outcome_cols

    # select outcome columns
    outcome_gwas <- outcome_gwas[, ..old_outcome_cols]

    print("selected outcome colnames")
    print(head(outcome_gwas))

    # change outcome colnames
    for (old_name in names(outcome_cols)) {

        if (outcome_cols[[old_name]] == "samplesize" & old_name == "") {
            outcome_gwas$samplesize <- opt$outcome_samplesize
        } else {
            setnames(outcome_gwas, old_name, outcome_cols[[old_name]])
        }
    }

    print("changed outcome colnames")
    print(head(outcome_gwas))

    # remove invalid entries
    outcome_gwas <- outcome_gwas[!is.na(outcome_gwas$pval),]

    print("valid outcomes")
    print(head(outcome_gwas))

    # limit each gwas to common variants
    outcome_gwas <- data.frame(outcome_gwas[outcome_gwas$SNP %in% exposure_gwas$SNP,])
    exposure_gwas <- exposure_gwas[exposure_gwas$SNP %in% outcome_gwas$SNP,]

    print("outcome after intersection")
    print(head(outcome_gwas))

    # use TwoSampleMR::format_data to get standardised data structures
    exposure_gwas <- format_data(exposure_gwas, type = "exposure", min_pval = 0)
    outcome_gwas <- format_data(outcome_gwas, type = "outcome", min_pval = 0)

    print("outcome after reformatting")
    print(head(outcome_gwas))

    print("exposure after reformatting")
    print(head(exposure_gwas))

    if (nrow(exposure_gwas) > 1) {
        clump <- plink_clump(
            exposure_gwas, 
            clump_kb = opt$plink_clump_kb, 
            clump_r2 = opt$plink_clump_r2, 
            clump_p = 1,
            linkage_file = opt$plink_linkage_files, 
            plink_bin = opt$plink_bin,
            plink_memory = opt$plink_memory
        )

        exposure_gwas <- exposure_gwas[which(exposure_gwas$SNP %in% clump$SNP),]
    }

    print("exposure gwas after plink clumping")
    print(head(exposure_gwas))

    outcome_gwas <- outcome_gwas[outcome_gwas$SNP %in% exposure_gwas$SNP,]

    print("outcome gwas after plink clumping")
    print(head(outcome_gwas))

    exposure_gwas$type <- "trans"

    if (opt$gene_positions != "") {
        
        gene_position_map <- data.frame(fread(opt$gene_positions))
        gene_position_map <- gene_position_map[which(gene_position_map$chromosome_name %in% 1:22),]
        act_map <- gene_position_map[tolower(gene_position_map$hgnc_symbol) == exposure_name,]
        
        # define which locus is a cis-pqtl if any - variant within gene/cis region
        if (nrow(act_map) > 0) {
            
            exposure_gwas$type[exposure_gwas$chr.exposure == act_map$chromosome_name &
                               exposure_gwas$pos.exposure > act_map$start_position - opt$cis_region & 
                               exposure_gwas$pos.exposure < act_map$end_position + opt$cis_region
            ] <- "cis"
        } else {
            message("gene not found in gene position file")
        }
        
    } else {
        message("gene_position_map_file is NULL; will find trans exposure variants only")
    }
    
    if (any(is.na(exposure_gwas$type))) {
        
        message("NA type for a clumped variant")
        return()  # error
    }

    all_harm <- exposure_gwas
    
    if (nrow(exposure_gwas) > 0) {
        all_harm$exposure <- paste(exposure_name, "all", sep="_")
        all_harm <- harmonise_data(all_harm, outcome_gwas, 1)
    } else {
        all_harm <- data.frame()
    }

    cis_harm <- exposure_gwas[exposure_gwas$type == "cis",]

    if (nrow(cis_harm) > 0) {
        cis_harm$exposure <- paste(exposure_name, "cis", sep="_")
        cis_harm <- harmonise_data(cis_harm, outcome_gwas, 1)
    } else {
        cis_harm <- data.frame()
    }

    harmonised_results <- list(cis_name=cis_harm, all_name=all_harm)
    names(harmonised_results) <- c(cis_name, all_name)

    return(harmonised_results)
}

files_out <- convert_input_gwas(opt)

# save exposure and outcome in harmonised tsmr format for cis, all
write.csv(files_out[[1]], paste0(names(files_out)[1], ".csv"), row.names = FALSE)
write.csv(files_out[[2]], paste0(names(files_out)[2], ".csv"), row.names = FALSE)
