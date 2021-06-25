// collate the results

process COLLATE_MR {

    label 'short'
    publishDir "$params.outdir/collated_mr/", mode: params.publish_dir_mode

    input:
    file mr_results
    file singlesnp
    file directionality
    file heterogeneity
    file pleiotropy

    output:
    path "*_mr_results.csv"
    path "*_singlesnp.csv"
    path "*_directionality.csv"
    path "*_heterogeneity.csv"
    path "*_pleiotropy.csv"
    path "*_full_results.csv"

    script:
    """
    echo $mr_results;
    echo $singlesnp;
    echo $directionality;
    echo $heterogeneity;
    echo $pleiotropy;
    
    Rscript  --verbose $baseDir/bin/collate_mr.R \
             --mr_results '$mr_results' \
             --singlesnp '$singlesnp' \
             --directionality '$directionality' \
             --heterogeneity '$heterogeneity' \
             --pleiotropy '$pleiotropy' \
             --auxiliary_script_dir '$baseDir/bin/auxiliary' \
             --num_exposures $params.num_exposures ;
    """
}

process COLLATE_COLOC {

    label 'short'
    publishDir "$params.outdir/collated_coloc/", mode: params.publish_dir_mode

    input:
    file coloc_res
    file coloc_susie_res

    output:
    path "*_coloc_res.csv"
    path "*_coloc_susie_res.csv"

    script:
    """
    echo $coloc_res;
    echo $coloc_susie_res;
    
    Rscript  --verbose $baseDir/bin/collate_coloc.R \
             --coloc_res '$coloc_res' \
             --coloc_susie_res '$coloc_susie_res' ;
    """
}
