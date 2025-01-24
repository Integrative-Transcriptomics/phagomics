#!/usr/bin/env nextflow

process generateOutputCsv {
    debug true
    publishDir "$params.outDir", mode: 'copy'

    input:
    path(resultUK)
    path(resultK)
}