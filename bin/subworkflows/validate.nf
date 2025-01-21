#!/usr/bin/env nextflow

include { foldseekValidate }        from "../modules/foldseek.nf"
include { generateFoldseekReport }  from "../modules/foldseek.nf"

workflow validateFoldseek {
    take:
    ch_structures

    main:
    ch_structures
    | foldseekValidate
    | generateFoldseekReport
    | collectFile( name: 'report.tsv', storeDir: "$params.outDir" )

    emit:
    'report.tsv'
}