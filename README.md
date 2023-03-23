# ONTmetacriptom-NF
A Nextflow Pipeline for Oxford Nanapore Technology(ONT) Long Reads Metatrascriptomic Data Analysis

[![nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.01.0-brightgreen.svg)](http://nextflow.io)

## Introduction 
**ONTmetacriptom-NF** is a bioinformatics analysis pipepliene used for the 
analysis of [Oxford Nanopore Technology](https://nanoporetech.com/) 
long reads metatrascriptomic data. The pipeline is built using [Nextflow](https://www.nextflow.io/), a workflow tool to run tasks across multiple
compute infrastructure in a very portable manner. It uses [docker](https://www.docker.com/) and [singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html) contantainers making installation inconsequential 
and results highly reproducible.

## Pipeline summary
**ONTmetacriptom-NF** pipeline performs by default the following:
- FastQ reads quality control
- Adaptor removal
- Sorting of ribosomal RNA and Messenger RNA
- Clustering of metatranscriptomes
- Quantification of transcripts and counting reads mapping to host

## Components 
**ONTmetacriptom-NF** uses the following software components and tools: 
* python=3.9.13
* pip=22.1.2
* NanoPlot 1.40.0
* mutiqc 1.13
* porechop 0.2.4
* samtools 1.6
* sortmerna 4.3.4
* isonclust 0.0.6.1
* spoa 4.0.7
* isoncorrect 0.0.8
* minimap2 2.17
* NanoCount 1.0.0post6

## Quick Start
1. Install [Nexflow](https://www.nextflow.io/docs/latest/getstarted.html#installation) (>=22.04.5)
2. Install any of [Docker](https://docs.docker.com/engine/install/) or [Singularity](https://singularity-tutorial.github.io/01-installation/)
3. Download the pipeline and test it on a minimal dataset in a single command:
```
nextflow run main.nf
```

4. To run the pipeline with your dataset, use the command:
```
nextflow run main.nf --readsDir '/path/to/dir_Containing_your_barcodes'
```

### To run a pipeline from a GitHub repository, 
1. You can run it by specifying the project name as shown below to test the pipeline:
```
nextflow run eanbit-rt/ONTmetacriptom-NF
```
2. To run the pipeline with your dataset.
```
nextflow run eanbit-rt/ONTmetacriptom-NF --readsDir '/path/to/dir_Containing_your_barcodes'
```
It automatically downloads it and store in the $HOME/.nextflow folder.

Use the command info to show the project information, e.g.:
```
nextflow info eanbit-rt/ONTmetacriptom-NF
```
## Credits
We thank the following people for their extensive assistance in the development of this pipeline (in alphabetical order):
- [Fedrick Kebaso](https://github.com/fredrickkebaso)
- [Samuel Oduor](https://github.com/samordil)
- [Stephen Kuria](https://github.com/sephoh)
- [Manase Aloo](https://github.com/manasealoo)
- [Mark Tefero Kivumbi](https://github.com/tefer0)

## References
These literatures shed more light on the background theories and concepts that were adapted to advance our methodology.
1. Shakya, M., Lo, C. C., & Chain, P. S. G. (2019). Advances and challenges in metatranscriptomic analysis. Frontiers in Genetics, 10(SEP), 904. https://doi.org/10.3389/FGENE.2019.00904/BIBTEX
2. Sahlin, K., Sipos, B., James, P. L., & Medvedev, P. (2021). Error correction enables use of Oxford Nanopore technology for reference-free    transcriptome analysis. Nature Communications, 12(1). https://doi.org/10.1038/S41467-020-20340-8
3. Joshua Batson, Gytis Dudas, Eric Haas-Stapleton, Amy L Kistler, Lucy M Li, Phoenix Logan, Kalani Ratnasiri, Hanna Retallack (2021) Single mosquito metatranscriptomics identifies vectors, emerging pathogens and reservoirs in one assay eLife 10:e68353 https://doi.org/10.7554/eLife.68353
