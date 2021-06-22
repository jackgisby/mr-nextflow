#' Temporarily install a package
#'
#' Temporarily install packages that aren't present on the cluster.
#'
#' @param package_name
#' Name of package to be installed.
#'
#' @param dependencies
#' If TRUE, install package dependences as well.
#'
#' @details
#' Add 'path' to .libPaths, and be sure that it is not
#' at the first position, otherwise any other package during
#' this session would be installed into 'path'.
#'
#' @references
#' https://stackoverflow.com/questions/14896941/install-an-r-package-temporarily-only-for-the-current-session
#'
#' @export

tmp.install.packages <- function(
  package_name,
  dependencies = TRUE,
  repos = "http://cran.us.r-project.org",
  github = FALSE,
  pkg_folder = NULL,
  ...
) {
  if (is.null(pkg_folder)) {
    path <- tempdir()
  } else {
    path <- pkg_folder
  }

  firstpath <- .libPaths()[1]

  .libPaths(c(firstpath, path))

  if (github) {
    withr::with_libpaths(new = path, devtools::install_github(package_name, dependencies = dependencies, ...))
  } else {
    install.packages(package_name, dependencies = dependencies, lib = path, repos = repos, ...)
  }
}

#' getting run settings from file
#'
#' Turns a csv containing settings for the run into a list
#'
#' @param settings_file_location
#' The path to the settings file for conversion to list
#'
#' @return
#' A list containing run settings

get_settings_from_file <- function(settings_file_location) {
  
  settings_file <- read.table(settings_file_location, sep = ",", quote = "\"")
  run_settings <- list()

  for (i in 1:nrow(settings_file)) {
    run_settings[[settings_file[i, 1]]] <- settings_file[i, 2]
  }

  run_settings[["clump_kb"]] <- as.numeric(run_settings[["clump_kb"]])
  run_settings[["clump_r2"]] <- as.numeric(run_settings[["clump_r2"]])
  run_settings[["outcome_samplesize"]] <- as.numeric(run_settings[["outcome_samplesize"]])
  run_settings[["p_threshold"]] <- as.numeric(run_settings[["p_threshold"]])
  run_settings[["cis_region"]] <- as.numeric(run_settings[["cis_region"]])
  run_settings[["num_proteins"]] <- as.numeric(run_settings[["num_proteins"]])
  
  if ("to" %in% names(run_settings)) {
    run_settings[["chr"]] <- as.numeric(run_settings[["chr"]])
    run_settings[["from"]] <- as.numeric(run_settings[["from"]])
    run_settings[["to"]] <- as.numeric(run_settings[["to"]])
    run_settings[["look_around"]] <- as.numeric(run_settings[["look_around"]])
    run_settings[["s"]] <- as.numeric(run_settings[["s"]])
    run_settings[["UKB_samplesize"]] <- as.numeric(run_settings[["UKB_samplesize"]])
  }
  
  if ("exposure_samplesize" %in% names(run_settings)) {
    run_settings[["exposure_samplesize"]] <- as.numeric(run_settings[["exposure_samplesize"]])
  }

  return(run_settings)
}

#' Get a list of BED files for use by plink
#'
#' Gets BED files or uses plink2 to convert VCFs to BEDs.
#'
#' @references
#' 1000G VCFs from FTP: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
#'
#' @export

get_bed_files <- function(linkage_file_dir, linkage_file_type = "UKB") {
  linkage_file_list <- list()

  for (f in list.files(linkage_file_dir)) {
    if (grepl(".bed", f)) {
      f <- gsub(".bed", "", f)

      if (linkage_file_type == "1000G") {
        chr <- gsub(".eur.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes", "", gsub("ALL.chr", "", f))
      } else if (linkage_file_type == "UKB") {
        chr <- gsub(".*_chr(\\d+)_10000_.*", "\\1", f)
      } else {
        chr <- gsub("_10000_random_unrelated_white_british", "", gsub("ukbb_chr", "", f))
      }
      
      if (as.integer(chr) %in% 1:22) {
        linkage_file_list[[as.integer(chr)]] <- paste0(linkage_file_dir, f)
      }
    }
  }
  
  return(linkage_file_list)
}
