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

**NanoPlot-** Plotting tool for long read sequencing data and alignments.Details can be accesed here [NanoPlot](https://github.com/wdecoster/NanoPlot).

**Multiqc-** MultiQC searches a given directory for analysis logs and compiles a HTML report. It's a general use tool, perfect for summarising the output from numerous bioinformatics tools.More information available here [multiqc](https://github.com/ewels/MultiQC).


