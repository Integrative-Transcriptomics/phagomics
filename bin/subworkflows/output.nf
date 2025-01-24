#!/usr/bin/env nextflow

include { clusterMembers }          from "../modules/cluster.nf"
include { generateFoldseekReport }  from "../modules/foldseek.nf"

workflow outputReport {
    take:
    ch_alnfile
    clufile

    main:
    clusterMembers( clufile )
    | first // convert to value-ch so it can be used indefinetly in generateFoldseekReport
    | set { clu_out }

    generateFoldseekReport( ch_alnfile, clu_out )
    | collectFile( name: 'report.tsv', storeDir: "$params.outDir", keepHeader: true )

    emit:
    'report.tsv'
}