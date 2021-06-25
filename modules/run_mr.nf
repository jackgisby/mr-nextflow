// run MR

process RUN_MR {

    label 'short'
    publishDir "$params.outdir/mr/", mode: params.publish_dir_mode

    input:
    file harm
    file mrinput

    output:
    path "*_mr_results.csv", emit: mr_results
    path "*_singlesnp.csv", emit: singlesnp
    path "*_directionality.csv", emit: directionality
    path "*_heterogeneity.csv", emit: heterogeneity
    path "*_pleiotropy.csv", emit: pleiotropy

    script:
    """
    echo $harm;
    echo $mrinput;

    Rscript  --verbose $baseDir/bin/run_mr.R \
             --harm '$harm' \
             --mrinput '$mrinput' \
             --auxiliary_script_dir '$baseDir/bin/auxiliary' ;
    """
}
