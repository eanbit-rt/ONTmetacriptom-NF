#!/usr/bin/env nextflow

/*
 * Copyright (c) 2022, EANBIT Residential training.
 * Copyright (c) 2022, Internation Centre of Insect and physiology (icipe).
 */
 
/* 
 * 'ONTmetacriptom-NF' - A Nextflow pipeline for ONT long reads 
 * metatranscriptomic data analysis.
 * 
 * This pipeline that reproduces steps from the GATK best practics of SNP 
 * calling with RNAseq data procedure:
 * https://software.broadinstitute.org/gatk/guide/article?id=3891
 * 
 * Fedrick Kebaso
 * Samuel Oduor
 * Stephen Kuria 
 */

/* 
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl = 2

/*
 * Define the default parameters
 */ 
params.reads = "$projectDir/data/*.gz"
params.outdir = "ONTresults"

log.info """\

    O N T m e t a c r i p t o m e - N F  v 0.1 
    ===========================================
    reads    : $params.reads
    results  : $params.outdir
"""
.stripIndent()

/* 
 * Import modules 
 
 include { 
    // Processes, functions and channels
  } from './modules/metatranscriptome.nf.nf' 
*/

/* 
 * main pipeline logic
 */

 workflow {
    // Section 1a: Quality Check using nanoqc


    /* Section 1b: Generating final report using MultiQC
     and outputs from 1a
    */


    // Section 2: ONT adaptor trimming using porechop tool


    // Section 3a: Downloading rRNA database
    

    /* Section 3b: Sorting mRNA from trRNA using ref rRNA from 
    sectiion 3a
    */
 }