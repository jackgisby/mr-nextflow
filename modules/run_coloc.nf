// selects variants (e.g. using plink clumping)

process RUN_COLOC {

    label 'long'
    publishDir "$params.outdir/coloc/", mode: params.publish_dir_mode

    input:
    file exposure_data
    file outcome_data
    file LD

    output:
    path "*_coloc_res.csv", emit: coloc_res
    path "*_coloc_obj.rds"
    path "*_coloc_susie_res.csv", emit: coloc_susie_res
    path "*_coloc_susie_obj.rds"

    script:
    """
    echo $exposure_data;
    echo $outcome_data;
    echo $LD;

    Rscript  --verbose $baseDir/bin/run_coloc.R \
             --exposure_input_data '$exposure_data' \
             --outcome_input_data '$outcome_data' \
             --LD '$LD' \
             --auxiliary_script_dir '$baseDir/bin/auxiliary' \
             --exposure_type '$params.coloc.exposure.type' \
             --exposure_s '$params.coloc.exposure.s' \
             --exposure_sdy '$params.coloc.exposure.sdy' \
             --outcome_type '$params.coloc.outcome.type' \
             --outcome_s '$params.coloc.outcome.s' \
             --outcome_sdy '$params.coloc.outcome.sdy' ;
    """
}
