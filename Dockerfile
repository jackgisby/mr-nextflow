# docker build -t jackgisby/mr-nextflow .
# docker push jackgisby/mr-nextflow

FROM r-base

RUN R -e "install.packages(c('data.table', 'optparse', 'TwoSampleMR', 'R.utils'))"
