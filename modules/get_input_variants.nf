// converts each exposure into TwoSampleMR dataset with the exposure
// only keeps significant variants to shrink file sizes

process GET_INPUT_VARIANTS {

    label 'plink_short'
    publishDir "$params.outdir/harmonised_data/", mode: params.publish_dir_mode

    input:
    file exposure_input_data

    output:
    path "*_cis_harm.csv", emit: cis_harm
    path "*_all_harm.csv",  emit: all_harm
    path "*_cis_exposure.csv", emit: cis_exposure
    path "*_cis_outcome.csv", emit: cis_outcome
    path "*_top_exposure.csv", emit: top_exposure
    path "*_top_outcome.csv", emit: top_outcome

    script:  // many arguments so we can pass the column names of these files
    """
    echo $exposure_input_data ;

    Rscript  --verbose $baseDir/bin/convert_input_gwas.R \
             --exposure_input_data '$exposure_input_data' \
             --outcome_input_data '$params.outcome.input_data' \
             --p_cutoff $params.p_cutoff \
             --auxiliary_script_dir '$baseDir/bin/auxiliary' \
             --cis_region '$params.cis_region' \
             --plink_memory $params.plink.clump_memory \
             --plink_clump_r2 $params.plink.clump_r2 \
             --plink_clump_kb $params.plink.clump_kb \
             --plink_bin '$params.plink.bin' \
             --plink_linkage_files '$params.plink.linkage_files' \
             --gene_positions '$params.gene_to_positions' \
             --gene_filenames '$params.gene_to_filenames' \
             --exposure_snp_col '$params.exposure.snp_col' \
             --exposure_beta_col '$params.exposure.beta_col' \
             --exposure_se_col '$params.exposure.se_col' \
             --exposure_eaf_col '$params.exposure.eaf_col' \
             --exposure_effect_allele_col '$params.exposure.effect_allele_col' \
             --exposure_other_allele_col '$params.exposure.other_allele_col' \
             --exposure_pval_col '$params.exposure.pval_col' \
             --exposure_samplesize_col '$params.exposure.samplesize_col' \
             --exposure_chr_col '$params.exposure.chr_col' \
             --exposure_pos_col '$params.exposure.pos_col' \
             --outcome_snp_col '$params.outcome.snp_col' \
             --outcome_beta_col '$params.outcome.beta_col' \
             --outcome_se_col '$params.outcome.se_col' \
             --outcome_eaf_col '$params.outcome.eaf_col' \
             --outcome_effect_allele_col '$params.outcome.effect_allele_col' \
             --outcome_other_allele_col '$params.outcome.other_allele_col' \
             --outcome_pval_col '$params.outcome.pval_col' \
             --outcome_samplesize_col '$params.outcome.samplesize_col' \
             --outcome_chr_col '$params.outcome.chr_col' \
             --outcome_pos_col '$params.outcome.pos_col' \
             --outcome_samplesize $params.outcome.samplesize ;
    """
}
