#!/usr/bin/env nextflow

include { foldseekReport }          from "../modules/scripts.nf"

workflow REPORT_NEW {
    take:
        alnfile

    main:
        foldseekReport( alnfile )
        | view
        // clusterReport()
        // postulatedReport()
        // interproscanReport()

}