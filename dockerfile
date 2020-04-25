FROM rocker/verse:3.6.3

ENV PKGS \
  openjdk-11-jdk\
  libxml2-dev \
  libcurl4

ENV WHEN 2020-04-25

RUN apt-get update
RUN apt-get install -y $PKGS

COPY /code home/rstudio/code
COPY /data home/rstudio/data

RUN R -e "options(repos = \
  list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/${WHEN}')); \
  install.packages(c('here', 'magrittr', 'pROC', 'hmeasure','pdp', 'rworldmaps',\
  'DALEX', 'ingredients', 'caret', 'e1071', 'icosa', 'ggthemes', \
  'patchwork', 'ggalluvial', 'raster', 'ggpubr', 'RCurl'))"

RUN R -e "if (! ('methods' %in% rownames(installed.packages()))) { install.packages('methods') };\
  if (! ('statmod' %in% rownames(installed.packages()))) { install.packages('statmod') }; \
  if (! ('stats' %in% rownames(installed.packages()))) { install.packages('stats') };\
  if (! ('graphics' %in% rownames(installed.packages()))) { install.packages('graphics') }; \
  if (! ('RCurl' %in% rownames(installed.packages()))) { install.packages('RCurl') }; \
  if (! ('jsonlite' %in% rownames(installed.packages()))) { install.packages('jsonlite') }; \
  if (! ('tools' %in% rownames(installed.packages()))) { install.packages('tools') }; \
  if (! ('utils' %in% rownames(installed.packages()))) { install.packages('utils') }"

RUN R -e "install.packages('h2o', type='source', \
  repos='http://h2o-release.s3.amazonaws.com/h2o/rel-tibshirani/2/R')"
