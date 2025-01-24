#!/usr/bin/env nextflow

include { foldseekValidate }        from "../modules/foldseek.nf"
include { generateValidationReport }  from "../modules/foldseek.nf"

workflow validateFoldseek {
    take:
    ch_structures

    main:
    ch_structures
    | foldseekValidate
    | generateValidationReport
    | collectFile( name: 'validationReport.tsv', storeDir: "$params.outDir", keepHeader: true )

    emit:
    'validationReport.tsv'
}