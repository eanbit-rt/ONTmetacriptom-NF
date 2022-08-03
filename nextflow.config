#!/opt/apps/nextflow/22.04.0/bin/nextflow

//Enable DSL 2 syntax

//Import modules here

include { nanoqc} from "./modules/qualitycheck.nf" addParams(outdir: "${params.outdir}")
include { porechop} from "./modules/qualitycheck.nf" addParams(outdir: "${params.outdir}")
include {concatenate} from "./modules/qualitycheck.nf" addParams(outdir: "${params.outdir}")
include {SORTMERNA} from "./modules/sortmeRNA.nf" addParams(outdir: "${params.outdir}")
include {MINIMAP2} from "./modules/minimap_tool.nf" addParams(outdir: "${params.outdir}")
// set the reads channel

Channel.fromPath( params.reads )
     .set{ read_ch }
Channel.fromPath(params.porechop)
    .set{porechop_ch}

Channel.fromPath(params.chopped)
    .set{choppedFiles_ch}

Channel.fromPath(params.refSeq)
    .set{refSeq_ch}
Channel.fromPath(params.sorted_reads)
    .set{sorted_reads_ch}




workflow{
    // Quality check and trimming
    // process 1a
    //nanoqc(read_ch)
    //SORTMERNA(concatenate(porechop(read_ch)))

    MINIMAP2(refSeq_ch,sorted_reads_ch)

    // porechop_ch.view()
    //concatenate(porechop_ch.collect())

    //SORTMERNA(choppedFiles_ch)

}