#!/usr/bin/env nextflow

include { foldseekReport }          from "../modules/scripts.nf"
include { clusterReport }           from "../modules/scripts.nf"
include { postulatedReport }        from "../modules/scripts.nf"


workflow REPORT_NEW {
    take:
        alnfile
        clusterfile
        proteinDescriptionsfile

    main:
        foldseekReport( alnfile )
        | view{it -> "Foldseek report: $it"}
        
        clusterReport( clusterfile, proteinDescriptionsfile )
        | flatten
        | map{ it -> [it.name[0..-9], it] } // restore [id, path] convention
        | view{it -> "Cluster report: $it"}
        
        postulatedReport( proteinDescriptionsfile )
        | flatten
        | map{ it -> [it.name[0..-9], it] } // restore [id, path] convention
        | take( 2 )
        | view{it -> "Postulated report: $it"}
        
        // interproscanReport()

}