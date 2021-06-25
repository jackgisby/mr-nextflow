# load libraries
library("optparse")
library("data.table")

# do argument parsing
option_list = list(
    make_option(c("--coloc_res"), type="character", default=NULL),
    make_option(c("--coloc_susie_res"), type="character", default=NULL)
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

print("options")
print(opt)

collate_coloc <- function(opt) {

    coloc_res <- data.frame(readr::read_csv(opt$coloc_res))
    print("coloc_res")
    print(head(coloc_res))

    coloc_susie_res <- data.frame(readr::read_csv(opt$coloc_susie_res))
    print("coloc_susie_res")
    print(head(coloc_susie_res))

    # if (ncol(coloc_susie_res) == 1) {
    #     coloc_susie_res <- data.frame()
    # }

    if (all(grepl("top", coloc_res$exposure))) {
        mr_run <- "top"
    } else if (all(grepl("cis", coloc_res$exposure))) {
        mr_run <- "cis"
    } else {
        stopifnot(FALSE)
    }

    return_names <- c(paste0(mr_run, "_coloc_res"), paste0(mr_run, "_coloc_susie_res"))

    final_results <- list("coloc_res"=coloc_res, "coloc_susie_res"=coloc_susie_res)
    names(final_results) <- return_names

    return(final_results)
}

files_out <- collate_coloc(opt)

for (i in 1:length(files_out)) {
    write.csv(files_out[[i]], paste0(names(files_out)[i], ".csv"), row.names = FALSE)
}
