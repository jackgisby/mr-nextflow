# load libraries
library("optparse")
library("data.table")

# load functions

# do argument parsing
option_list = list(
    make_option(c("--exposure_input_data"), type="character", default=NULL, 
                help="exposure input file name", metavar="character"),
    make_option(c("--outcome_input_data"), type="character", default=NULL, 
                help="outcome input file name", metavar="character"),
    make_option(c("--data_format_exposure"), type="character", default=NULL, 
                help="the file extension for the exposure", metavar="character"),
    make_option(c("--data_format_outcome"), type="character", default=NULL, 
                help="the file extension for the outcome", metavar="character"),
    make_option(c("--p_cutoff"), type="numeric", default=NULL, 
                help="p value cutoff for the exposure", metavar="character")
)
 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# load exposure w/ fread
exposure_gwas <- fread(opt$exposure_input_data)
exposure_gwas <- exposure_gwas[exposure_gwas$RSID != ".",]
exposure_gwas <- exposure_gwas[exposure_gwas$RSID != "",]
exposure_gwas <- exposure_gwas[!is.na(exposure_gwas$PVAL),]
print(exposure_gwas)

# limit exposure to significant variants
exposure_gwas <- data.frame(exposure_gwas[exposure_gwas$PVAL <= opt$p_cutoff,])
# exposure_gwas$CHR_POS <- paste(exposure_gwas$CHR, exposure_gwas$POS, sep="_")

print(exposure_gwas)

# load outcome w/ fread
outcome_gwas <- fread(opt$outcome_input_data)
outcome_gwas <- outcome_gwas[!is.na(outcome_gwas$Pval),]
outcome_gwas <- outcome_gwas[outcome_gwas$RSID != "",]
print(outcome_gwas)

# limit each gwas to common variants
# outcome_gwas$CHR_POS <- paste(exposure_gwas$CHR, exposure_gwas$POS, sep="_")
outcome_gwas <- data.frame(outcome[outcome_gwas$RSID %in% exposure_gwas$RSID,])
exposure_gwas <- data.frame(outcome[exposure_gwas$RSID %in% outcome_gwas$RSID,])

# use TwoSampleMR::format_data to get standardised data structures
outcome_gwas <- format_data(
    outcome_gwas, type = "outcome", 
    chr_col = "CHR", pos_col = "POS", snp_col = "RSID", beta_col = "Effect", 
    se_col = "SE", eaf_col = "AF Effect Allele", effect_allele_col = "Effect Allele",
    other_allele_col = "Non-Effect Allele", pval_col = "Pval"
)

exposure_gwas <- format_data(
    exposure_gwas, type = "exposure", 
    snp_col = "RSID", beta_col = "BETA", se_col = "SE", eaf_col = "FREQ1",
    effect_allele_col = "EFFECT_ALLELE", other_allele_col = "REFERENCE_ALLELE",
    pval_col = "PVAL", samplesize_col = "N", chr_col = "CHR", pos_col = "POS"
)

# save exposure and outcome in tsmr format
write.csv(exposure_gwas, "exposure_gwas.csv")
write.csv(outcome_gwas, "outcome_gwas.csv")
