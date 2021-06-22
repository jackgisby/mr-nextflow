# load libraries
library("optparse")
library("data.table")
library("TwoSampleMR")

# do argument parsing
option_list = list(
    make_option(c("--exposure_input_data"), type="character", default=NULL, help="exposure input file name", metavar="character"),
    make_option(c("--outcome_input_data"), type="character", default=NULL, help="outcome input file name", metavar="character"),
    make_option(c("--data_format_exposure"), type="character", default=NULL, help="the file extension for the exposure", metavar="character"),
    make_option(c("--data_format_outcome"), type="character", default=NULL, help="the file extension for the outcome", metavar="character"),
    make_option(c("--auxiliary_script_dir"), type="character", default=NULL, help="the location of helper scripts", metavar="character"),
    make_option(c("--p_cutoff"), type="numeric", default=NULL),
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

# save exposure and outcome in tsmr format
write.csv(exposure_gwas, "exposure_gwas.csv")
write.csv(outcome_gwas, "outcome_gwas.csv")
