#!/usr/bin/env nextflow

process generateOutputCsv {
    debug true
    publishDir "$params.outDir", mode: 'copy'
    containerOptions { "--rm" }

    input:
    path(resultUK)
    path(resultK)
}