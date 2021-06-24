// selects variants (e.g. using plink clumping)

process RUN_COLOC {

    label 'plink_long'
    publishDir "$params.outdir/coloc/", mode: params.publish_dir_mode

    input:
    file exposure_data
    file outcome_data

    output:
    path "*_coloc_res.csv"
    path "*_coloc_susie_abf.csv"
    path "*_LD.csv"

    script:
    """
    echo $exposure_data;
    echo $outcome_data;
    
    Rscript  --verbose $baseDir/bin/run_coloc.R \
             --exposure_input_data '$exposure_data' \
             --outcome_input_data '$outcome_data' \
             --auxiliary_script_dir '$baseDir/bin/auxiliary' \
             --exposure_type '$params.coloc.exposure.type' \
             --exposure_s '$params.coloc.exposure.s' \
             --exposure_sdy '$params.coloc.exposure.sdy' \
             --outcome_type '$params.coloc.outcome.type' \
             --outcome_s '$params.coloc.outcome.s' \
             --outcome_sdy '$params.coloc.outcome.sdy' \
             --plink_bin '$params.plink.bin' \
             --plink_linkage_files '$params.plink.linkage_files' \
             --plink_memory $params.plink.matrix_memory \
             --nref $params.coloc.nref ;
    """
}
