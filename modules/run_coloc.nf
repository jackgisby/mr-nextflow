// selects variants (e.g. using plink clumping)

process RUN_COLOC {

    label 'long_plink'

    input:
    file exposure_input_data
    file location

    output:
    path "*_coloc_res.csv"
    path "*_coloc_susie_res.csv"
    
    script:
    """
    echo $exposure_input_data;
    echo $location;
    
    Rscript  --verbose $baseDir/bin/run_coloc.R \
             --exposure_input_data '$exposure_input_data' \
             --genomic_region '$location' ;
    """
}
