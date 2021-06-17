// converts each exposure into TwoSampleMR dataset with the exposure
// only keeps significant variants to shrink file sizes

process GET_INPUT_VARIANTS {

    label 'get_gwas_data'

    input:
    file exposure_input_data

    output:
    file "exposure_gwas.tsv", emit: tsmr_exposure
    file "outcome_gwas.tsv",  emit: tsmr_outcome
    
    script:  // must use argument parser in R and must also allow input of colnames
    """
    R convert_input_gwas.R 
        --exposure_input_data $exposure_input_data 
        --outcome_input_data $params.outcome_input_data 
        --data_format_exposure $params.data_format.exposure 
        --data_format_outcome $params.data_format.outcome 
        --p_cutoff $params.p_cutoff
    """
}
