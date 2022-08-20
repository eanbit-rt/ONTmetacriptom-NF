# --------------------------------------------------------------------
# Creates a docker image with all the necessary tools and software   -
# needed by the ONTmetacriptome-NF pipeline.			     -
# --------------------------------------------------------------------

FROM mambaorg/micromamba
LABEL Your name='Samuel Oduor'

RUN \
  micromamba install -y -n base -c defaults -c bioconda -c conda-forge \
    isonclust=0.0.6.1 \
    isoncorrect=0.0.8 \
    minimap2=2.17 \
    nanoqc=0.9.4 \
    porechop=0.2.4 \
    pycoqc=2.5.2 \
    sortmerna=4.3.4 \
    nanocount==1.0.0 \
  && micromamba clean -a -y
