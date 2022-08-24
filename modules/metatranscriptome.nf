#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Process 1A: Checking the quality of the ONT reads 
 * with nanqc
 */
process NANOPLOT_QC {
    publishDir "${params.outdir}/qcOutput", mode: 'copy'
    tag "Quality Check"
    
    input:
      tuple(val(fastqDir),  path(fastqFile))
    output:
        path "${fastqDir}"

    script:
    """
    NanoPlot \
    --threads 2 \
    --outdir ${fastqDir} \
    --prefix ${fastqDir} \
    --fastq ${fastqFile}
    """
}

 /*
 * Process 1b:  Collecting the outputs of nanoqc process 
 * to create a final report using MultiQC tool.
 */
 process MULTIQC_REPORT {
    publishDir "${params.outdir}/multiqcOutput", mode:'copy'
    tag 'QC Report Aggregate'

    input:
    path '*'
    
    output:
    path 'aggregated_QC_report.html'
    
    script:
    """
    multiqc --filename aggregated_QC_report.html . 
    """
}

/*
* Process 2: ONT adapter removal using porechop tool
* finds and removing adapters from Oxford Nanopore reads.
*/
process  PORECHOP_TRIM {
    publishDir "${params.outdir}", mode: 'copy'
    tag "Trimming Reads"

  input:
    path readsDir
  
  output:
    path 'Poreched_Concat_Dir'
  
  script:
  """
  porechop \
  --input ${readsDir} \
  --format fastq \
  --barcode_dir Poreched_Concat_Dir
  """
}


/*
* Process 3a: Downloading reference rRNA database from
* https://github.com/biocore/sortmerna/archive/2.1b.zip
* using download_rnadb.sh a custome script in bin 
* directory
*/


/*
* Process 3b: filter rRNA fragments from metatranscriptomic data
* using the downloaded rRNA database
*/

/*
* Process 4: clustering long transcriptomic reads into gene 
* families using isOnCLust tool
*/

/*
* Process 5: error-correcting Oxford Nanopore cDNA reads. 
* It is designed to handle highly variable coverage and 
* exon variation within reads using isOnCorrect tool
*/

/*
*Process 6: 
* minimap2
* sequence alignment program that aligns DNA or mRNA 
* sequences against a large reference database
*/

/*
* Process 7:
* Estimates transcript abundances using NanoCount
*/

/*
* Process 8
* HTSeq is a Python package that calculates the number 
* of mapped reads to each gene.
*/
