#!/usr/bin/env nextflow     

Channel.fromPath(params.reads)
    .view()


