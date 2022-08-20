#!/usr/bin/env bash

# Install pycoqc for checking ont read quality
conda install -c aleg -c anaconda -c bioconda -c conda-forge pycoqc=2.5.2

# Install nanoqc for nanopore quality check
#Investigate nucleotide composition and base quality.
conda install -c bioconda nanoqc

# Install porechop for trimming ont adapters
conda install -c bioconda porechop

# Install sortmerna
conda install -c bioconda sortmerna

# Install isonclust tool for 
conda install -c bioconda isonclust

# Install isoncorrect 
conda install -c bioconda isoncorrect

# Install 
conda install -c bioconda minimap2

# Install Nanocout tool for abundance
# First installation
pip install git+https://github.com/a-slide/NanoCount.git

pip install git+https://github.com/a-slide/NanoCount.git --upgrade
