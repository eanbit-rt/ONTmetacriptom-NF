# Nextflow Pipeline for ONT Long Reads Metatrascriptomic Data Analysis

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

## References
1. These literatures shed more light on the background theories and concepts that were adapted to advance our methodology.
Shakya, M., Lo, C. C., & Chain, P. S. G. (2019). Advances and challenges in metatranscriptomic analysis. Frontiers in Genetics, 10(SEP), 904. https://doi.org/10.3389/FGENE.2019.00904/BIBTEX
2. Sahlin, K., Sipos, B., James, P. L., & Medvedev, P. (2021). Error correction enables use of Oxford Nanopore technology for reference-free    transcriptome analysis. Nature Communications, 12(1). https://doi.org/10.1038/S41467-020-20340-8
3. Joshua Batson, Gytis Dudas, Eric Haas-Stapleton, Amy L Kistler, Lucy M Li, Phoenix Logan, Kalani Ratnasiri, Hanna Retallack (2021) Single mosquito metatranscriptomics identifies vectors, emerging pathogens and reservoirs in one assay eLife 10:e68353 https://doi.org/10.7554/eLife.68353
