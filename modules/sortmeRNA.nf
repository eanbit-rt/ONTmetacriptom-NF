#!/usr/bin/env/ nextflow

// Running DSL2
nextflow.enable.dsl = 2

//Params.outdir = "${baseDir}/sortmernaOutput"

// choppedFiles_ch = Channel.fromPath(
//     "${baseDir}/porechop_output/*.fastq", checkIfExists: true)


process SORTMERNA {
    publishDir "${params.outdir}/sortmernaOutput", mode:'copy'

    input:
        file(chopped)

    output:
        path "rRNA_sequences"
        path "mRNA_sequences"

    script:
        filename = chopped.simpleName
        // out = "rRNABlastOutput"
        // log = "rRNALogOutput"

        """
        mkdir -p rRNA_sequences mRNA_sequences

        sortmerna \
        --ref ${PWD}/rRNA_databases/silva-arc-16s-id95.fasta \
        --workdir ${filename} \
        --reads ${chopped} \
        --aligned rRNA_sequences/${filename}_sorted_rRNA \
        --fastx --other mRNA_sequences/${filename}_sorted_mRNA \
        --blast '1 cigar qcov' ${filename} 

        
        """
}

// workflow {
//     result_ch = SORTMERNA(choppedFiles_ch)
// }