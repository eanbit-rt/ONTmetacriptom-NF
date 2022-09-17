#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Process 1A: Checking the quality of the ONT reads 
 * with NanoPlot
 * Takes a list of sequences per barcode and processes 
 * them to produce single reports/stats per barcode
 */
process NANOPLOT_QC {
    tag "Pre-trim Quality Check"
    
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
    publishDir "${params.outdir}", mode:'copy'
    tag 'Pre-trim QC Report Aggregate'

    input:
    path '*'
    
    output:
    path 'Pre_Trim_QC_report.html'
    
    script:
    """
    multiqc --filename Pre_Trim_QC_report.html . 
    """
}

/*
* Process 2: ONT adapter removal using porechop tool
* finds and removing adapters from Oxford Nanopore reads.
*/
process  PORECHOP_TRIM {
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

// Post Adaptor Removal Quality Check
process POST_TRIM_NANOPLOT_QC {
    tag "Post-trim Quality Check"
    
    input:
      path trimmedFastqFile
    output:
        path "${dirName}_without_Adapters"

    script:
      dirName = trimmedFastqFile.simpleName
    """
    NanoPlot \
    --threads 4 \
    --outdir "${dirName}_without_Adapters" \
    --prefix ${dirName} \
    --fastq ${trimmedFastqFile}
    """
}

// Post Adopter Removal MultiQC
process POST_TRIM_MULTIQC_REPORT {
    publishDir "${params.outdir}", mode:'copy'
    tag 'QC Report Aggregate'

    input:
    path '*'
    
    output:
    path 'Post_Trim_QC_report.html'
    
    script:
    """
    multiqc --filename Post_Trim_QC_report.html . 
    """
}

/*
* Process 3a: Downloading reference rRNA database from
* https://github.com/biocore/sortmerna/archive/2.1b.zip
* using download_rnadb.sh a custom script in bin 
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
    --ref ${rnaDB[0]} \
    --ref ${rnaDB[1]} \
    --ref ${rnaDB[2]} \
    --ref ${rnaDB[3]} \
    --ref ${rnaDB[4]} \
    --ref ${rnaDB[5]} \
    -reads ${trimmedReadFile} \
    --workdir $PWD/work \
    --kvdb kvdb \
    -threads 4 \
    --aligned seqrRNA/${filename} \
    --fastx --other seqmRNA/${filename} \
    
  """
  }

/*
* Process 4: clustering long transcriptomic reads into gene 
* families using isOnCLust tool
*/
process ISONCLUST {
  tag "Gene Families Clustering"

  input:
    path seqReads
    
  output:
  path "${geneFamilies}"

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
process ISONCORRECT {
  tag "Reads Error correction"

  input:
    path(readsDir)
  
  output:
    path "${name}.fastq"

  script:
    name = readsDir.simpleName

  """
  run_isoncorrect \
  --fastq_folder ${readsDir} \
  --t 3 \
  --k 13 \
  --w 20 \
  --outfolder "${name}"

  cat ${name}/*/corrected_reads.fastq \
  > "${name}.fastq"
"""
}

/*
*Process 6: 
* minimap2
* sequence alignment program that aligns DNA or mRNA 
* sequences against a large reference database
*/
// 6a: Downloading the reference genome sequence for minimap2
process REFSEQ_GCF_016920715_DOWNLOAD {
    tag "Reference Seq Genome Downlooad"

    output:
        path 'GCF_016920715.1*.idx'

    script:
        """
        wget \
        https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/920/715/GCF_016920715.1_AaraD3/GCF_016920715.1_AaraD3_genomic.fna.gz
        
        gunzip GCF_016920715.1_AaraD3_genomic.fna.gz

        minimap2 \
        -d GCF_016920715.1_AaraD3_genomic.idx GCF_016920715.1_*
        """
}
// 6b: Mapping to the reference genome sequence
process MINIMAP2 {
  tag "Mapping to ref genome"

    input:
    path correctedReads
    path(refSeqIndexed)

    output:
    path "*.bam"

    script:
    filename = correctedReads.simpleName

    """
    minimap2 \
    -t 4 \
    -ax map-ont \
    -p 0 \
    -N 10 ${refSeqIndexed} ${correctedReads} | \
    samtools view -bh | \
    samtools sort > "${filename}.bam"
    """
}

/*
* Process 7:
* Estimates transcript abundances using NanoCount
*/
process NANOCOUNT {
  tag "transcripts count"
      publishDir "${params.outdir}", mode:'copy'

  input:
    path bam_file

  output:
    path "*.tsv"

  script:
    filename = bam_file.simpleName

  """
  NanoCount \
  --alignment_file ${bam_file} \
  --count_file "${filename}.trascripts.tsv" \
  --sec_scoring_threshold 0.8 \
  --verbose
  """
}

/*
* Process 8
* HTSeq is a Python package that calculates the number 
* of mapped reads to each gene.
* Ongoing process
*/