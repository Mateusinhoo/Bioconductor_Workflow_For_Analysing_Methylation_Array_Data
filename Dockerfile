# Dockerfile for methyCleanr
FROM rocker/r-bioc:3.18

RUN R -q -e "install.packages(c('yaml','ggplot2','matrixStats','R.utils'), repos='https://cloud.r-project.org')"
RUN R -q -e "BiocManager::install(c('minfi','limma','DMRcate','missMethyl','IlluminaHumanMethylation450kmanifest','IlluminaHumanMethylation450kanno.ilmn12.hg19','IlluminaHumanMethylationEPICmanifest','IlluminaHumanMethylationEPICanno.ilm10b4.hg19'), ask=FALSE, update=FALSE)"

WORKDIR /work
COPY . /work

CMD ["Rscript","-e","source('run_all.R')"]
