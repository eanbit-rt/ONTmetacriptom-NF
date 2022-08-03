#!/usr/bin/env/ nextflow

nextflow.enable.dsl = 2

//params.outdir = "${baseDir}/minimapOuput"

// refSeq_ch = Channel.fromPath("${baseDir}/minmapRefSeq/*.fa")
// reads_ch = Channel.fromPath("${baseDir}/all_corrected_reads/*.fastq")

process MINIMAP2 {
    publishDir "${params.outdir}/minimapOuput", mode:'copy'

    input:
    each sorted_reads
    path refSeq

    output:
    path "*.sam"

    script:
    filename = sorted_reads.simpleName

    """
    
}
// workflow{
//     MINIMAP2(reads_ch, refSeq_ch)
// }