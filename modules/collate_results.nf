// collate the results

process COLLATE_MR {

    label 'short'

    input:
    file harm

    output:
    path "exposure_variants_clumped.csv", emit: tsmr_exposure_clumped
    
    script:
    """
    echo $exposure_variants;
    
    Rscript  --verbose $baseDir/bin/collate_results.R \
             --exposure_input_data '$exposure_variants' ;
    """
}
