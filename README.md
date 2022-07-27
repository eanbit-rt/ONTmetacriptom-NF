# Nextflow Pipeline for ONT Long Reads Metatrascriptomic data analysis

## Introduction 
**ontmetacriptom.nf** is a bioinformatics analysis pipepliene used for the 
analysis of [Oxford Nanopore Technology](https://nanoporetech.com/) 
long reads metatrascriptomic data. The pipeline is built using [Nextflow](https://www.nextflow.io/), a workflow tool to run tasks across multiple
compute infrastructure in a very portable manner. It uses [docker](https://www.docker.com/) and [singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html) contantainers making installation inconsequential 
and results highly reproducible.

## Pipeline summary
**ontmetacriptom.nf** pipeline performs by default the following:
- Sequencing quality control
- Trimming of reads


## Quick Start
1. Install [Nexflow](https://www.nextflow.io/docs/latest/getstarted.html#installation) (>=22.04.5)
2. Install any of [Docker](https://docs.docker.com/engine/install/) or [Singularity](https://singularity-tutorial.github.io/01-installation/)
3. Download the pipeline and test it on a minimal dataset in a single command:
```
nextflow run ontmetacriptom.nf
```
## Credits
We thank the following people for their extensive assistance in the development of this pipeline (in alphabetical order):
- [Fedrick Kebaso](https://github.com/fredrickkebaso)
- [Manase Aloo](https://github.com/manasealoo)
- [Mark Kivumbi](https://github.com/tefer0)
- [Samuel Oduor](https://github.com/samordil)
- [Stephen Kuria](https://github.com/sephoh)
