# mr-nextflow

MR for multiple exposures implemented using nextflow. Plink clumping is used for instrument selection using the exposure. Various parameters can be changed, which can be found in `nextflow.config` - all of these can be changed from default and some are required to be changed from default, see `tests/test.config`. Specific columns are expected for both the outcome and the exposure, however the names of these columns can be specified in the config file. 

By default, multiple variant MR is carried out using both the TwoSampleMR package and the MendelianRandomization package (which takes into account variant correlations using the LD matrix). Colocalisation testing is done for the strongest (smallest p value) exposure and the cis region of the exposure (if a gene or a protein). 

Plink is used extensively, so the location of plink (v1.9) and bfile must be given.

If the pipeline is to be run on a cluster, a config file must be given and singularity used. This will likely need to overwrite some of the options in `conf/base.config`.

This workflow is, at present, a skeleton to be developed further. Currently, it will likely require adjustments for new data. 
