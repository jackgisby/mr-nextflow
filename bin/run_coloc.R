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

get_blank_return <- function() {

    susie_coloc_result_colnames <- c("nsnps", "hit1", "hit2", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf", "idx1", "idx2", "exposure", "to", "from", "region_width")
    susie_coloc_summary <- data.frame(matrix(ncol = length(susie_coloc_result_colnames), nrow=0))
    colnames(susie_coloc_summary) <- susie_coloc_result_colnames

    coloc_result_colnames <- c("PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf", "p12", "pass", "exposure", "to", "from", "region_width")
    coloc_summary <- data.frame(matrix(ncol = length(coloc_result_colnames), nrow=0))
    colnames(coloc_summary) <- coloc_result_colnames

    return(list("coloc_res"=coloc_summary, "coloc_obj"=data.frame(), "coloc_susie_res"=susie_coloc_summary, "coloc_susie_obj"=data.frame()))
}

run_coloc <- function(opt) {

    region_name <- gsub(".csv", "", gsub("_exposure", "", opt$exposure_input_data))

    exposure <- data.frame(data.table::fread(opt$exposure_input_data))
    print("exposure")
    print(head(exposure))

    outcome <- data.frame(data.table::fread(opt$outcome_input_data))
    print("outcome")
    print(head(outcome))

    LD <- as.matrix(data.table::fread(opt$LD))
    print("LD")
    print(LD)

    null_return <- get_blank_return()
    names(null_return) <- paste(region_name, names(null_return), sep="_")

    if (nrow(outcome) < 2 | nrow(exposure) < 2) {
        return(null_return)
    }

    rownames(LD) <- colnames(LD)

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
    if (typeof(susie_coloc_res) == "list" & length(susie_coloc_res) == 3) {
        susie_coloc_res <- process_sensitivity_data(susie_coloc_res$summary, region_name, to, from, nrow(exposure_in))
    } else {
        susie_coloc_res <- get_blank_return()[["coloc_susie_res"]]
    }

    # set up final results
    final_results <- list("coloc_res"=standard_sensitivity, "coloc_obj"=standard_coloc, "coloc_susie_res"=susie_coloc_res, "coloc_susie_obj"=susie_coloc_res)
    names(final_results) <- paste(region_name, names(final_results), sep="_")

    return(final_results)
}

files_out <- run_coloc(opt)

# save exposure and outcome in harmonised tsmr format for cis, all

for (i in 1:length(files_out)) {
    if (grepl("obj", names(files_out)[i])) {
        saveRDS(files_out[[i]], paste0(names(files_out)[i], ".rds"))
    } else {
        write.csv(files_out[[i]], paste0(names(files_out)[i], ".csv"), row.names = FALSE)
    }
}
