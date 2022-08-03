#!/usr/bin/env/ nextflow

// Running DSL2
//params.outdir = "${baseDir}/sortmernaOutput"

//choppedFiles_ch = Channel.fromPath(
    //"${baseDir}/porechop_output/*.fastq", checkIfExists: true)


process SORTMERNA {
    //publishDir params.outdir, mode:'copy'
    publishDir "${params.outdir}/sortmernaOutput/", mode:'copy'

    input:
        //file(pore_chopped_file)
        file(chopped)

    output:
        path "rRNABlastOutput"
        path "rRNALogOutput"

    script:
        //filename = pore_chopped_file.simpleName
        filename = chopped.simpleName

        // out = "rRNABlastOutput"
        // log = "rRNALogOutput"

        """
        mkdir -p rRNABlastOutput rRNALogOutput

        sortmerna \
        --ref ${PWD}/rRNA_databases/silva-arc-16s-id95.fasta \
        --workdir ${PWD}/${filename} \
        --reads ${chopped} \
        --aligned rRNA_sequences/${filename}_sorted_rRNA \
        --fastx --other mRNA_sequences/${filename}_sorted_mRNA \
        --blast '1 cigar qcov' ${filename} 

        mv rRNA_sequences/*.blast rRNABlastOutput
        mv rRNA_sequences/*.log rRNALogOutput

        rm -rf "${PWD}/${filename}/sortmerna/run/kvdb"
        """
}

// workflow {
//      result_ch = SORTMERNA(choppedFiles_ch)
//  }
