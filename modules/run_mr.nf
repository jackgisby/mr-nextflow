// run MR

process RUN_MR {

    label 'short'

    input:
    file harm
    file mr_input
    file matrix

    output:
    path "*_mr_results.csv", emit: mr_results
    path "*_top_location.csv", emit: top_location
    
    script:
    """
    echo $exposure_variants;
    
    Rscript  --verbose $baseDir/bin/run_mr.R \
             --exposure_input_data '$exposure_variants' ;
    """
}
