#!/usr/bin/env nextflow

include { interproscan } from "../modules/interproscan.nf" 

workflow INTERPROSCAN {
    take:
        unknownProteins // fasta file
    
    main:
        unknownProteins
        | splitFasta( by: 30 )
        | take( 15 )
        | interproscan

    // emit:
    //     results
}