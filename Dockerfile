FROM ubuntu:16.04

ADD package.json /home/package.json
ADD requirements.txt /home/requirements.txt
RUN apt-get update && apt-get install -y \
    vim \
    htop \
    unzip \
    curl \
    build-essential \
    libcurl4-openssl-dev \
    git \
    tmux \
    libxml2-dev \
    python3-setuptools \
    libssl-dev && \
    cd /home/ && \
    curl -sL https://deb.nodesource.com/setup_6.x | bash - && \
    apt-get install -y \
    nodejs \
    r-base \
    r-base-dev \
    python3-pip && \
    pip3 install -r /home/requirements.txt && \
    echo "install.packages('devtools', repos='http://cran.rstudio.com/')" | Rscript - && \
    echo "devtools::install_version('readr', version = '1.0.0', repos = 'http://cran.us.r-project.org')" | Rscript - && \
    echo "devtools::install_version('dplyr', version = '0.5.0', repos = 'http://cran.us.r-project.org')" | Rscript - && \
    echo "devtools::install_version('argparse', version = '1.0.4', repos = 'http://cran.us.r-project.org')" | Rscript - && \
    echo "devtools::install_version('purrr', version = '0.2.2', repos = 'http://cran.us.r-project.org')" | Rscript - && \
    echo "devtools::install_version('stringr', version = '1.2.0', repos = 'http://cran.us.r-project.org')" | Rscript - && \
    echo "devtools::install_version('lazyeval', version = '0.2.0', repos = 'http://cran.us.r-project.org')" | Rscript - && \
    echo "source('https://bioconductor.org/biocLite.R'); biocLite('GSEABase')" | Rscript - && \
    npm install && \
    apt-get remove -y \
    build-essential && \
    apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /home
