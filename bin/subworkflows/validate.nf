#!/usr/bin/env nextflow

include { foldseekValidate }            from "../modules/foldseek.nf"
include { generateValidationReport }    from "../modules/scripts.nf"

workflow VALIDATE {
    take:
        knownProteins

    main:
        knownProteins
        | foldseekValidate
        | generateValidationReport
        | collectFile( name: 'validationReport.tsv', storeDir: "$params.outDir", keepHeader: true )
}