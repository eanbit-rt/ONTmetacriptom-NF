/*
 * Process 1A: Checking the quality of the ONT reads 
 * with nanqc
 */
 

 /*
 * Process 1b:  Collecting the outputs of nanoqc process 
 * to create a final report using MultiQC tool.
 */


/*
* Process 2: ONT adaptor removal using porechop tool
* finds and removing adapters from Oxford Nanopore reads.
*/


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
* Process 8:
* HTSeq is a Python package that calculates the number 
* of mapped reads to each gene.
*/