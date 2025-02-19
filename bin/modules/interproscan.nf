#!/usr/bin/env nextflow

process interproscan {
    debug true
    publishDir "$params.outDir", mode: 'copy'
    
    input:
    val(unknownProteins)


    script:
    """
    """
}