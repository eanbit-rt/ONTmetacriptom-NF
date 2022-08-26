#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def rnadb_available( path ) {
/*
Searches for 'rRNA_database' dir
in the project directory and returns
true when found
*/
  myDir = file(path)
  allFiles = myDir.list()

    for( def file : allFiles ) {
        if( file == 'rRNA_databases') {
            return true
        }
        else {
            return false
            }
    }
}

/*
 * Process 1A: Checking the quality of the ONT reads 
 * with nanqc
 */
process NANOPLOT_QC {
    publishDir "${params.outdir}/qcOutput", mode: 'copy'
    tag "Quality Check"
    
    input:
      tuple(val(name),  path(fastqFile))
    output:
        path "${name}"

    script:
    """
    NanoPlot \
    --threads 4 \
    --outdir ${name} \
    --prefix ${name} \
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
    path 'General_QC_report.html'
    
    script:
    """
    multiqc --filename General_QC_report.html . 
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
    path 'porechopOuput'
  
  script:
  """
  porechop \
  --input ${readsDir} \
  --format fastq \
  --barcode_dir porechopOuput \
  --discard_unassigned
  """
}

/*
* Process 3a: Downloading reference rRNA database from
* https://github.com/biocore/sortmerna/archive/2.1b.zip
* using download_rnadb.sh a custome script in bin 
* directory
*/
rnaDatabase = "${projectDir}"
process GET_RNADATABASE {
    tag "Searching for rRNA database "
    publishDir "${projectDir}}", mode: 'copy', overwrite: true
    
    output:
        path 'rRNA_databases'

    script:
    if (rnadb_available(rnaDatabase))
        """
        mkdir rRNA_databases
        """
    else
        """
        getRNAdb.sh
        """
}

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
