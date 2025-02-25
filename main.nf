#!/usr/bin/env nextflow

include { FILTER                } from "./bin/subworkflows/filter.nf"
include { CLUSTER               } from "./bin/subworkflows/cluster.nf"
include { STRUCTURE_PREDICITON  } from "./bin/subworkflows/structure.nf"
include { SEARCH                } from "./bin/subworkflows/search.nf"
include { INTERPROSCAN          } from "./bin/subworkflows/interproscan.nf"
include { VALIDATE              } from "./bin/subworkflows/validate.nf"
include { REPORT                } from "./bin/subworkflows/output.nf"
include { REPORT_NEW            } from "./bin/subworkflows/output_new.nf"


workflow {

    ///
    /// IMPORT & FILTERING
    ///

    // Import protein files and limit length to 1500 residues
    Channel.fromPath( params.proteins, checkIfExists: true )
    | splitFasta( record: [id:true, desc:true, seqString:true] )
    | filter { record -> record.seqString.length() < 100 }
    | map{ it -> [id:it.id.replace("lcl|", ""), desc:it.desc, seqString:it.seqString]} // remove lcl| from start of id
    | set { ch_allProteins }

    FILTER( ch_allProteins )


    ///
    /// INPERPROSCAN
    ///

    INTERPROSCAN( FILTER.out.unknownProteins )


    ///
    /// CLUSTERING
    ///

    CLUSTER( FILTER.out.allProteins )


    ///
    /// STRUCTURE PREDICTION
    ///

    STRUCTURE_PREDICITON( CLUSTER.out.splitClusterReps )


    ///
    /// FOLDSEEK SEARCH
    /// 

    // Branch known/unknown
    STRUCTURE_PREDICITON.out
    | branch { it ->
        known:      it =~ /.*_known_.*/
        unknown:    it =~ /.*_unknown_.*/
    }
    | set { structures }

    SEARCH( structures.unknown )


    ///
    /// VALIDATION & OUTPUT
    /// 

    VALIDATE( structures.known )
    REPORT( SEARCH.out , CLUSTER.out.clusterMembers )

    REPORT_NEW( 
        SEARCH.out, 
        CLUSTER.out.clusterMembers,
        FILTER.out.proteinDescriptions
    )
}   

workflow.onComplete {
    println ( workflow.success ? """
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        outDir      : ${params.outDir}
        exit status : ${workflow.exitStatus}
        """ : """
        Failed:       ${workflow.errorReport}
        exit status : ${workflow.exitStatus}
        """
    )
}