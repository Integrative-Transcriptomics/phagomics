#!/usr/bin/env nextflow

include { foldseek }      from "../modules/foldseek.nf"

workflow SEARCH {
    take:
        unknownProteins
    
    main:
        // Run foldseek with unknown proteins 
        unknownProteins
        | foldseek
        | set{ alignments }
    
    emit:
        alignments
}