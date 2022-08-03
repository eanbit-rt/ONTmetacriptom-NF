//Process 1a: Quality check --Tool: Nanoqc


process nanoqc{
    publishDir "${params.outdir}/nanoqc_raw", mode:"copy", overwrite:false

    tag "Quality checking raw reads"

    input:

    path(reads)

    output:
    
    path(bar)

    script:
    
    bar = reads.baseName

      
    """
    
    nanoQC -o ${bar} ${reads}

    """
}
//porechop -i ${reads} --format fastq -o ${bar}
//cat "./Results/trimmed/${bar}"  >> porechopped.fastq
//file "porechopped.fastq"
//    porechop -i ${reads} --format fastq -o ${bar}
//cat ${PWD}/${params.outdir}/trimmed/*.fastq  >> porechopped.fastq
    

process porechop{

    publishDir "${params.outdir}/trimmed", mode:"copy", overwrite:false
    
    tag "Trimming raw reads"

    input:
    path(reads)

    output:
    path(bar)
    
    script:

    bar=reads.baseName

    """
    porechop -i ${reads} --format fastq -o ${bar}


    """
}

process concatenate{

    publishDir "${params.outdir}", mode:"copy"
    
    tag "Concatenating porechopped sequences"

    input:
    path porechop

    output:
    path "porechopDir/porechopped.fastq"

    script:
        
    """
    
    mkdir  porechopDir

    echo ${porechop}
    cat ${porechop} >> porechopDir/porechopped.fastq
    
    """
}


