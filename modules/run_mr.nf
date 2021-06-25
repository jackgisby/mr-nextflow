// run MR

process RUN_MR {

    label 'short'

    input:
    file harm
    file mr_input

    output:
    path "*_mr_results.csv", emit: mr_results
    
    script:
    """
    echo $exposure_variants;
    
    Rscript  --verbose $baseDir/bin/run_mr.R \
             --harm '$harm' \
             --mr_input '$mr_input';
    """
}
