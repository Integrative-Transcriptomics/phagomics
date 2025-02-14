#!/usr/bin/env nextflow

include { foldseek }      from "../modules/foldseek.nf"

workflow SEARCH {
    take:
        unknownProteins
    
    main:
        // Run unknown proteins with foldseek
        unknownProteins
        | foldseek
        | set{ alignments }
    
    emit:
        alignments
}