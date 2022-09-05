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
params.readsDir = "$projectDir/data"
params.outdir = "ONTresults"
params.refGenome = ''

log.info """\

    O N T m e t a c r i p t o m e - N F  v 0.1 
    ===========================================
    readsDir    : $params.readsDir
    results     : $params.outdir
"""
.stripIndent()

/* 
 * Import modules 
*/

 include { 
      NANOPLOT_QC;
      MULTIQC_REPORT;
      PORECHOP_TRIM;
      DOWNLOAD_rRNADATABASE;
      SORTMERNA;
      ISONCLUST;
      ISONCORRECT;
      REFSEQ_GCF_016920715_DOWNLOAD;
      MINIMAP2;
      NANOCOUNT
  } from './modules/metatranscriptome.nf' 


/* 
 * main pipeline logic
 */
 workflow {
    channel
      .fromPath("${params.readsDir}/**/*.gz")
      .map { fastq -> tuple(fastq.parent.name, fastq)}
      .groupTuple()
      .set { raw_reads_ch }

/*
 the channel is converted into value channell to enable 
 multiple use of the same values without the channel 
 consuming it
*/
    channel
    .value(file("${projectDir}/rRNA_databases/*.fasta"))
    .collect()
    .set { rRNAdatabases }

/*
    channel
    .value(file("${params.refGenome}"))
    .set { refSeq_ch }
*/
    // Section 1a: Quality Checking
    NANOPLOT_QC(raw_reads_ch)
    
    /* Section 1b: Generating final report using 
    * outputs from 1a
    */
    MULTIQC_REPORT(NANOPLOT_QC.out.collect())
  
    // Section 2: ONT adaptor removal
    PORECHOP_TRIM(params.readsDir)
    //PORECHOP_TRIM.out.flatten().view()

    // Section 3a: Download rRNA database. Need internet connection
    //DOWNLOAD_rRNADATABASE()

    // Section 3b: rRNA fragments filtering
    SORTMERNA(//DOWNLOAD_rRNADATABASE.out.collect(),
              rRNAdatabases,
              PORECHOP_TRIM.out.flatten()
              )

    /* 
    *Section 4: Clustering genes into falmilies using 
    *isONclust
    */
    ISONCLUST(SORTMERNA.out)

    //Seciton 5: ONT reads error correction
    ISONCORRECT(ISONCLUST.out
    .map{ dirpath -> tuple(dirpath.name, dirpath) }
    )

    // Section 6a: Download a reference datatabase
    REFSEQ_GCF_016920715_DOWNLOAD()

    // Section 6: Alignment to a reference datatabase
    MINIMAP2(ISONCORRECT.out, 
            REFSEQ_GCF_016920715_DOWNLOAD.out
            )
    // Section 7: Transcript abundance estimation
    NANOCOUNT(MINIMAP2.out)
    //NANOCOUNT.out.view()

    // section 8: Calculating the number of mapped reads to each gene
 }