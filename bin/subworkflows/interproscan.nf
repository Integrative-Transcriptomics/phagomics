#!/usr/bin/env nextflow

include { interproscan } from "../modules/interproscan.nf" 

workflow INTERPROSCAN {
    take:
        unknownProteins
    
    main:
        interproscan( unknownProteins )

    // emit:
    //     results
}