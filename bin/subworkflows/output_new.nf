#!/usr/bin/env nextflow

include { foldseekReport }          from "../modules/scripts.nf"
include { clusterReport }           from "../modules/scripts.nf"
include { postulatedReport }        from "../modules/scripts.nf"
include { mergeReports }            from "../modules/scripts.nf"


workflow REPORT_NEW {
    take:
        alnfile
        clusterfile
        proteinDescriptionsfile
        interproscanReport

    main:
        foldseekReport( alnfile )
        | map{ it -> tuple((it.name =~ /^(.*?)_prot/)[0][1] , it) } // [phage, prot_id, path]
        | set{ fs }
        
        clusterReport( clusterfile, proteinDescriptionsfile )
        | flatten
        | map{ it -> tuple((it.name[0..-9] =~ /^(.*?)_prot/)[0][1] , it) } // [phage, prot_id, path]
        | set{ cl }
        
        postulatedReport( proteinDescriptionsfile )
        | flatten
        | map{ it -> tuple((it.name[0..-9] =~ /^(.*?)_prot/)[0][1] , it) } // [phage, prot_id, path]
        | set{ ps }
        
        // interproscanReport is done in modules/interproscan.nf
        interproscanReport
        | flatten
        | map{ it -> tuple((it.name[0..-9] =~ /^(.*?)_prot/)[0][1] , it) } // [phage, prot_id, path]
        | set{ ip }

        fs.concat ( cl, ps, ip )
        | groupTuple()
        | mergeReports

}