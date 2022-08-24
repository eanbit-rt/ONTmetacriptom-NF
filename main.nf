#!/usr/bin/env nextflow

/*
 * Copyright (c) 2022, EANBIT Residential training.
 * Copyright (c) 2022, International Centre of Insect Physiology 
 * and Ecology (icipe).
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
params.reads = "$projectDir/data/barcode*/*.gz"
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
*/

 include { 
      QUALITY_NANOPLOT;
      //REPORT_MULTIQC
  } from './modules/metatranscriptome.nf' 


/* 
 * main pipeline logic
 */

 workflow {
    channel
      .fromPath(params.reads)
      .map { fastq -> tuple(fastq.parent.name, fastq)}
      .groupTuple()
      .set { raw_reads_ch }


    // Section 1a: Quality Checking
    QUALITY_NANOPLOT(raw_reads_ch)
    
    /* Section 1b: Generating final report using 
    * outputs from 1a
    */
   //REPORT_MULTIQC(QUALITY_NANOQC.out.collect())
   //REPORT_MULTIQC.out.view()


    // Section 2: ONT adaptor removal


    // Section 3a: Downloading rRNA database
    

    // Section 3b: rRNA fragments filtering
 
    // Section 4: Clustering genes into falmilies

    // Seciton 5: ONT reads error correction

    // Section 6: Alignment to a reference datatabase

    // Section 7: Transcript abundance estimation

    // section 8: Calculating the number of mapped reads to each gene
 }