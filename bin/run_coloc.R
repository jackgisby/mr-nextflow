# load libraries
library("optparse")

# do argument parsing
option_list = list(
    make_option(c("--exposure_input_data"), type="character", default=NULL),
    make_option(c("--location"), type="character", default=NULL)
)
 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

print("options")
print(opt)

# load exposure w/ fread
exposure_gwas <- read.csv(opt$exposure_input_data)

# do plink clumping
exposure_gwas <- exposure_gwas[sample(nrow(exposure_gwas), 5),]

# write clumped output data
write.csv(exposure_gwas, "exposure_variants_clumped.csv")
