#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Process 1A: Checking the quality of the ONT reads 
 * with NanoPlot
 * Takes a list of sequences per barcode and processes 
 * then to produce single reports/stats per barcode
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
    tag "Trim Reads"

  input:
    path readsDir
  
  output:
    path 'porechopOuput/*.fastq'
  
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
process DOWNLOAD_rRNADATABASE {
    tag "Download rRNA ref databases"

    output:
        path 'rRNA_databases/*'

    script:
        """
        getRNAdb.sh
        """
}

/*
* Process 3b: filter rRNA fragments from metatranscriptomic data
* using the downloaded rRNA database
*/
process SORTMERNA {
    publishDir "${params.outdir}/sortmernaOutput", mode:'copy'
    tag 'Sorting'
  input:
    path rnaDB
    path trimmedReadFile
  
  output:
    path "seqmRNA/*"
  
  script:
    filename = trimmedReadFile.simpleName

  """
    sortmerna \
    --workdir $PWD/work \
    --kvdb $filename/kvdb \
    --ref ${rnaDB[0]} \
    --ref ${rnaDB[1]} \
    --ref ${rnaDB[2]} \
    --ref ${rnaDB[3]} \
    --ref ${rnaDB[4]} \
    --ref ${rnaDB[5]} \
    -reads ${trimmedReadFile} \
    --aligned seqrRNA/${filename}_rRNA \
    --fastx --other seqmRNA/${filename}_mRNA 
  """
  }

/*
* Process 4: clustering long transcriptomic reads into gene 
* families using isOnCLust tool
*/
process ISONCLUST {
  publishDir "${params.outdir}/isonclustOuput", mode:'copy'
  tag "Gene Families Clustering"

  input:
    path seqReads
    
  output:
    // path "isonclustOutput/*" // Output for the first command
  path(geneFamilies)

  script:
    geneFamilies = seqReads.simpleName

  """ 
  isONclust \
  --fastq ${seqReads} \
  --ont --outfolder "isonclustOutput/${geneFamilies}"

  isONclust write_fastq \
  --clusters "isonclustOutput/${geneFamilies}/final_clusters.tsv" \
  --fastq ${seqReads} \
  --outfolder ${geneFamilies} \
  --N 1
  """
}

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

// ============== START OF ALTERNATIVE CODES ============
/*
def rnadb_available( path ) {
/*
Searches for 'rRNA_database' dir
in the project directory and returns
true when found. Needs the activation of GET_RNADATABASE
process
*/ ///////////////////////////////////////

/*
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
*\

/*
* Process 3a
* Checks whether there is already an existing rRNA database
* in the projectDir

rnaDatabase = "${projectDir}"
process  GET_RNADATABASE {
    tag "Searching for rRNA database "
    publishDir "${projectDir}", mode: 'copy', overwrite: false
    
    output:
        path 'rRNA_databases/*'

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
*/
// ======================= END =========================