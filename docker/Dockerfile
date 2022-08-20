# --------------------------------------------------------------------
# Creates a docker image with all the necessary tools and software   -
# needed by the ONTmetacriptome-NF pipeline.			     -
# --------------------------------------------------------------------

FROM mambaorg/micromamba
LABEL name='Group 2 EANBIT icipe 2022'

RUN \
  micromamba install -y -n base -c conda-forge -c anaconda -c bioconda -c defaults \ 
    git pip \
    multiqc=1.10.1 \
    salmon=1.5.1 \
    fastqc=0.11.9 \
    pycoqc=2.5.2 \
    nanoqc=0.9.4 \
    porechop=0.2.4 \
    sortmerna=4.3.4 \ 
    isonclust=0.0.6.1 \
    isoncorrect=0.0.8 \
    minimap2=2.17 \ 
  && micromamba clean -a -y

RUN pip install git+https://github.com/a-slide/NanoCount.git
RUN pip install git+https://github.com/a-slide/NanoCount.git --upgrade