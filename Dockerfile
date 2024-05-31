FROM rocker/tidyverse:4.4.0

LABEL maintainer = "Ryan Corbett (corbettr@chop.edu)"


#########################################
RUN apt-get update && apt-get install -y --no-install-recommends apt-utils dialog

# Add curl, bzip2 and some dev libs
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    curl \
    bzip2 \
    zlib1g \
    libbz2-dev \
    liblzma-dev \
    libreadline-dev

# libmagick++-dev is needed for coloblindr to install
RUN apt-get -y --no-install-recommends install \
    libgdal-dev \
    libudunits2-dev \
    libmagick++-dev

# Set the Bioconductor repository as the primary repository
RUN R -e "options(repos = BiocManager::repositories())"

# Install BiocManager and the desired version of Bioconductor
RUN R -e "install.packages('BiocManager', dependencies=TRUE)"
RUN R -e "BiocManager::install(version = '3.19')"

# Install packages
RUN R -e 'BiocManager::install(c( \
  "circlize", \
  "ComplexHeatmap", \
  "data.table", \
  "ggthemes", \
  "readxl", \
  "RColorBrewer", \
  "survival", \
  "survMisc", \
  "survminer" \
))'
  
	
RUN R -e "remotes::install_github('clauswilke/colorblindr', ref = '1ac3d4d62dad047b68bb66c06cee927a4517d678', dependencies = TRUE)"
RUN R -e "remotes::install_github('thomasp85/patchwork', ref = '1cb732b129ed6a65774796dc1f618558c7498b66')"
	
WORKDIR /rocker-build/

ADD Dockerfile . 
