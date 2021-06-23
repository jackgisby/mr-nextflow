// uses plink to get LD matrix (and related harmonised dataset - flipped w/ reference to matrix)

process GET_LD_MATRIX {

    label 'short_plink'

    input:
    file harm

    output:
    path "*_mr_input.csv", emit: mr_input
    path "*_matrix.csv", emit: matrix
    
    script:
    """
    echo $exposure_variants;
    
    Rscript  --verbose $baseDir/bin/get_ld_matrix.R \
             --tsmr_harmonised '$harm' ;
    """
}
