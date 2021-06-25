# load libraries
library("optparse")
library("coloc")
library("susieR")

# do argument parsing
option_list = list(
    make_option(c("--exposure_input_data"), type="character", default=NULL, help="exposure input file name", metavar="character"),
    make_option(c("--outcome_input_data"), type="character", default=NULL, help="outcome input file name", metavar="character"),
    make_option(c("--LD"), type="character", default=NULL, help="linkage disequilibrium matrix", metavar="character"),
    make_option(c("--auxiliary_script_dir"), type="character", default=NULL, help="the location of helper scripts", metavar="character"),
    make_option(c("--exposure_type"), type="character", default=NULL),
    make_option(c("--exposure_s"), type="character", default=NULL),
    make_option(c("--exposure_sdy"), type="character", default=NULL),
    make_option(c("--outcome_type"), type="character", default=NULL),
    make_option(c("--outcome_s"), type="character", default=NULL),
    make_option(c("--outcome_sdy"), type="character", default=NULL)
)
 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

print("options")
print(opt)

source(paste0(opt$auxiliary_script_dir, "/coloc.R"))
source(paste0(opt$auxiliary_script_dir, "/linkage.R"))

run_coloc <- function(opt) {

    region_name <- gsub(".csv", "", gsub("_exposure", "", opt$exposure_input_data))

    exposure <- read.csv(opt$exposure_input_data)
    print("exposure")
    print(head(exposure))

    outcome <- read.csv(opt$outcome_input_data)
    print("outcome")
    print(head(outcome))

    LD <- as.matrix(read.csv(opt$LD))
    print("LD")
    print(LD)
    rownames(LD) <- colnames(LD)

    return_names <- c(paste0(region_name, "_coloc_res"), paste0(region_name, "_coloc_obj"), paste0(region_name, "_coloc_susie_abf"), paste0(region_name, "_coloc_susie_obj"))
    null_return <- list("standard_coloc"=data.frame(), "standard_obj"=data.frame(), "susie_abf"=data.frame(), "susie_obj"=data.frame())
    names(null_return) <- return_names

    if (nrow(outcome) == 0 & nrow(exposure) == 0) {
        return(null_return)
    }

    chr <- exposure$chr[1]
    to <- max(exposure$pos)
    from <- min(exposure$pos)
    
    # coloc format
    exposure_in <- get_coloc_list(exposure, list_type = opt$exposure_type, s = opt$exposure_s, sdy = opt$exposure_sdy)
    outcome_in <- get_coloc_list(outcome, list_type = opt$outcome_type, s = opt$outcome_s, sdy = opt$outcome_sdy)

    # do coloc and store sensitivity analysis
    standard_coloc <- coloc.abf(exposure_in, outcome_in)
    standard_sensitivity <- sensitivity(standard_coloc, rule = "H4 > 0.5")
    standard_sensitivity <- process_sensitivity_data(standard_sensitivity, region_name, to, from, nrow(exposure_in))

    exposure_in_ld <- get_coloc_list(exposure, list_type = opt$exposure_type, s = opt$exposure_s, sdy = opt$exposure_sdy, LD=LD)
    outcome_in_ld <- get_coloc_list(outcome, list_type = opt$outcome_type, s = opt$outcome_s, sdy = opt$outcome_sdy, LD=LD)

    exposure_run_susie <- runsusie(exposure_in_ld)
    outcome_run_susie <- runsusie(outcome_in_ld)

    susie_coloc_res <- coloc.susie(exposure_run_susie, outcome_run_susie)
    susie_coloc_abf <- susie_coloc_res$results

    # set up final results
    final_results <- list("standard_coloc"=standard_sensitivity, "standard_obj"=standard_coloc, "susie_abf"=susie_coloc_abf, "susie_obj"=susie_coloc_res)
    names(final_results) <- return_names

    return(final_results)
}

files_out <- run_coloc(opt)

# save exposure and outcome in harmonised tsmr format for cis, all

for (i in 1:length(files_out)) {
    if (i %in% c(2, 4)) {
        saveRDS(files_out[[i]], paste0(names(files_out)[i], ".rds"))
    } else {
        write.csv(files_out[[i]], paste0(names(files_out)[i], ".csv"), row.names = FALSE)
    }
}
