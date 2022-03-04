# docker build -t jackgisby/mr-nextflow .
# docker push jackgisby/mr-nextflow

# use the TSMR container as a base
FROM mrcieu/twosamplemr

# install extra packages
RUN R -e "install.packages(c('data.table', 'optparse', 'R.utils', 'devtools', 'MendelianRandomization', 'readr'))"
RUN R -e "devtools::install_github(c('chr1swallace/coloc'))"

# get plink
RUN apt-get update --allow-releaseinfo-change && apt-get install -y curl unzip && apt-get clean
RUN curl 'https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20210606.zip' > /tmp/plink.zip
RUN cd /tmp && unzip plink.zip && rm plink.zip
RUN mv /tmp/plink /usr/local/bin
