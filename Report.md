# EANBiT 2022 Residential Training Mini Project Report
## Project: Nextflow Pipeline for ONT Long Reads Metatrascriptomic data (ONTmetacriptom-NF).
### Introduction

High throughput sequencing (shotgun and ONT) has accelerated characterization of microbes and their functionality in environment. 
Most of tools and pipelines are optimized for short reads metatranscriptomics. There are few standardized workflow for metatranscriptomic analysis of long read. 
In this project we set out to develop a pipeline for Oxford Nanopore reads particularly for metatranscriptomics data.

### Objectives

- To develop a standardised workflow for Oxford Nanopore reads using nextflow workflow language.
- Expand the workflow to make use of containers and scale to the cloud.

### Tools incorporated in the workflow.

**Nanoqc -**  Quality control tools for long read sequencing data aiming to replicate some of the plots made by fastQC. For more information check [nanoqc](https://github.com/wdecoster/nanoQC).

**NanoPlot -** Plotting tool for long read sequencing data and alignments.Details can be accesed here [NanoPlot](https://github.com/wdecoster/NanoPlot).

**Multiqc -** MultiQC searches a given directory for analysis logs and compiles a HTML report. It's a general use tool, perfect for summarising the output from numerous bioinformatics tools.More information available here [multiqc](https://github.com/ewels/MultiQC).

**Porechop -** Porechop is a tool for finding and removing adapters from Oxford Nanopore reads. It's described in details here [porechop](https://github.com/rrwick/Porechop).

**SortMeRNA -** SortMeRNA is a local sequence alignment tool for filtering, mapping and clustering. Details are given here [sortmerna](https://github.com/biocore/sortmerna). 

**isONclust -** isONclust is a tool for clustering either PacBio Iso-Seq reads, or Oxford Nanopore reads into clusters, where each cluster represents all reads that came from a gene. Output is a tsv file with each read assigned to a cluster-ID. Detailed information is available here [isONclust](https://github.com/ksahlin/isONclust).

**isONcorrect -** isONcorrect is a tool for error-correcting Oxford Nanopore cDNA reads. It is designed to handle highly variable coverage and exon variation within reads and achieves about a 0.5-1% median error rate after correction. It leverages regions shared between reads from different isoforms achieve low error rates even for low abundant transcripts. Here is more information [isONcorrect](https://github.com/ksahlin/isONcorrect).

**Minimap2 -** Minimap2 is a versatile sequence alignment program that aligns DNA or mRNA sequences against a large reference database. Its described here [minimap2](https://github.com/lh3/minimap2).


