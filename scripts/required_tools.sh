#!/usr/bin/env bash

# Create a conda env
conda create -n mini-proj python=3 pip

# Install pycoqc for checking ont read quality
#conda install -c aleg -c anaconda -c bioconda -c conda-forge pycoqc=2.5.2

# Install nanoqc for nanopore quality check
#Investigate nucleotide composition and base quality.
#conda install -c bioconda nanoqc 
# Install NanoPlot for checking the read quality of ONT reads
pip install NanoPlot
pip install NanoPlot --upgrade # to upgrade to the latest version

# Install Multiqc for generating one report from the NanoPlot output
pip install multiqc

# Install porechop for trimming ont adapters
conda install -c bioconda porechop

# Install sortmerna
conda install -c bioconda sortmerna

# Install isonclust tool for 
conda install -c bioconda isonclust

# Install isoncorrect 
#Install isONcorrect and its dependency `spoa`.
pip install isONcorrect
conda install -c bioconda spoa

# Install 
conda install -c bioconda minimap2

# Install Nanocout tool for abundance
# First installation
pip install git+https://github.com/a-slide/NanoCount.git

pip install git+https://github.com/a-slide/NanoCount.git --upgrade