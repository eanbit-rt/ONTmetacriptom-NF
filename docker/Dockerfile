# -----------------------------------------------------------------------
# Creates a docker image with all the necessary tools needed to run ONT -
# metatrancriptomic data analysis.				        -
# -----------------------------------------------------------------------

# Base image, a micromamba environment.
FROM mambaorg/micromamba

LABEL version="v1.0.0"
LABEL name='Group 2 EANBIT icipe 2022'
LABEL author="samordil@gmail.com"

# Setting the working directory
WORKDIR /home/mambauser

# installing using micromamba
RUN micromamba install -y -n base \
 -c conda-forge \
 -c bioconda \
 -c defaults \
    python=3.9.13 \
    porechop=0.2.4 \
    samtools=1.6 \
    sortmerna=4.3.4 \
    isonclust=0.0.6.1 \
    minimap2=2.17 \
    spoa=4.0.7 \
    && micromamba clean -a -y

# Installing using pip
RUN /opt/conda/bin/pip install \
 NanoPlot==1.40.0 \
 multiqc==1.12 \ 
 isONcorrect==0.0.8 \
 NanoCount==1.0.0post6

# Clean after pip 
RUN /opt/conda/bin/pip cache purge
