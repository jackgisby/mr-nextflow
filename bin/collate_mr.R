# load libraries
library("optparse")
library("data.table")

# do argument parsing
option_list = list(
    make_option(c("--mr_results"), type="character", default=NULL),
    make_option(c("--singlesnp"), type="character", default=NULL),
    make_option(c("--directionality"), type="character", default=NULL),
    make_option(c("--heterogeneity"), type="character", default=NULL),
    make_option(c("--pleiotropy"), type="character", default=NULL),
    make_option(c("--mr_results_leaveoneout"), type="character", default=NULL),
    make_option(c("--directionality_leaveoneout"), type="character", default=NULL),
    make_option(c("--heterogeneity_leaveoneout"), type="character", default=NULL),
    make_option(c("--pleiotropy_leaveoneout"), type="character", default=NULL),
    make_option(c("--auxiliary_script_dir"), type="character", default=NULL, help="the location of helper scripts", metavar="character"),
    make_option(c("--gene_filenames"), type="integer", default=NULL)
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

print("options")
print(opt)

source(paste0(opt$auxiliary_script_dir, "/mr.R"))

collate_mr_results <- function(opt) {

    mr_results <- data.frame(readr::read_csv(opt$mr_results))
    print("mr_results")
    print(head(mr_results))

    singlesnp <- data.frame(readr::read_csv(opt$singlesnp))
    print("singlesnp")
    print(head(singlesnp))

    directionality <- data.frame(readr::read_csv(opt$directionality))
    print("directionality")
    print(head(directionality))

    heterogeneity <- data.frame(readr::read_csv(opt$heterogeneity))
    print("heterogeneity")
    print(head(heterogeneity))

    pleiotropy <- data.frame(readr::read_csv(opt$pleiotropy))
    print("pleiotropy")
    print(head(pleiotropy))

    mr_results_leaveoneout <- data.frame(readr::read_csv(opt$mr_results_leaveoneout))
    print("mr_results_leaveoneout")
    print(head(mr_results_leaveoneout))

    null_return <- list("_mr_results"=data.frame(), "_directionality"=data.frame(), "_heterogeneity"=data.frame(), "_pleiotropy"=data.frame(), "_mr_results_leaveoneout"=data.frame(), "_singlesnp"=data.frame(), "_full_results"=data.frame())

    if (nrow(mr_results) == 0) {
        return(null_return)
    }

    if (all(grepl("all", mr_results$exposure))) {
        mr_run <- "all"
    } else if (all(grepl("cis", mr_results$exposure))) {
        mr_run <- "cis"
    } else {
        stopifnot(FALSE)
    }
    
    full_results <- get_collated_df(mr_results, opt$gene_filenames, directionality, heterogeneity, pleiotropy)
    print("full results")
    print(head(full_results))

    singlesnp <- singlesnp[!grepl("ALL", singlesnp$SNP),]

    final_results <- list("mr_results"=mr_results, "directionality"=directionality, "heterogeneity"=heterogeneity, "pleiotropy"=pleiotropy, 
                          "mr_results_leaveoneout"=mr_results_leaveoneout, "singlesnp"=singlesnp, "full_results"=full_results)
    names(final_results) <- paste(mr_run, names(final_results), sep="_")
    return(final_results)
}

files_out <- collate_mr_results(opt)

for (i in 1:length(files_out)) {
    write.csv(files_out[[i]], paste0(names(files_out)[i], ".csv"), row.names = FALSE)
}
