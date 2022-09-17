# -----------------------------------------------------------------------
# Creates a docker image with all the necessary tools needed to run ONT -
# metatrancriptomic data analysis.				                        -
# -----------------------------------------------------------------------

# Base image
FROM debian:jessie

# File Author / Maintainer
LABEL maintainer="Samuel Oduor samordil@gmail.com"
LABEL description="This is custom Docker Image for \
    ONT long reads metatrascriptomic data analysis. \
    and their dependancies"

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH=/opt/conda/bin:$PATH

################## BEGIN INSTALLATION ######################
RUN apt-get update && apt-get install -y \
    wget \
    unzip \
    time \
&& \
    apt-get clean \
&& \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* \
&& \
    wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \ 
        -O /miniconda.sh \
&& \
    /bin/bash /miniconda.sh -b -p /opt/conda \
&& \
    rm /miniconda.sh \
&& \
    conda install -y -n base \
        -c conda-forge \
        -c bioconda \
        -c defaults \
        pip=22.1.2 \
        python=3.9.13 \
        porechop=0.2.4 \
        samtools=1.6 \
        sortmerna=4.3.4 \
        isonclust=0.0.6.1 \
        minimap2=2.17 \
        spoa=4.0.7 \
&& \ 
    conda clean -a -y \
&& \
    /opt/conda/bin/pip install \
        NanoPlot==1.40.0 \
        multiqc==1.13 \ 
        isONcorrect==0.0.8 \
        NanoCount==1.0.0post6 \
&& \
    /opt/conda/bin/pip cache purge \
&& \
    useradd -m docker && echo "docker:docker" | chpasswd && adduser docker sudo

USER docker

RUN echo "source activate" >> ~/.bashrc 

# To add a new user to the system
# RUN useradd --create-home --shell /bin/bash jessie && \
#   chown -R jessie:jessie /home/jessie
# USER jessie