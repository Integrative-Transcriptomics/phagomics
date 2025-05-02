#!/usr/bin/env nextflow

include { interproscan } from "../modules/interproscan.nf" 

workflow INTERPROSCAN {
    take:
        unknownProteins // fasta file
    
    main:
        unknownProteins
        | splitFasta( record: [id:true, seqString: true] )  // complicated way to filter out known proteins
        | filter{ record -> record.id =~ /_unknown$/ }      // because cluster reps are not filtered
        | map { record -> ">${record.id}\n${record.seqString}\n" }
        | collectFile
        | splitFasta( by:30 )
        | interproscan
        | set{ report }

    emit:
        report
}