# docker build -t jackgisby/mr-nextflow .
# docker push jackgisby/mr-nextflow

FROM mrcieu/twosamplemr

RUN R -e "install.packages(c('data.table', 'optparse', 'R.utils'))"
