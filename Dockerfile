# docker build -t jackgisby/mr-nextflow .
# docker push jackgisby/mr-nextflow

FROM mrcieu/twosamplemr

RUN R -e "install.packages(c('data.table', 'optparse', 'R.utils', 'devtools', 'MendelianRandomization', 'readr'))"
RUN R -e "devtools::install_github(c('chr1swallace/coloc'))"
