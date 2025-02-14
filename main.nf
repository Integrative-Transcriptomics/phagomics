#!/usr/bin/env nextflow

include { FILTER                } from "./bin/subworkflows/filter.nf"
include { CLUSTER               } from "./bin/subworkflows/cluster.nf"
include { STRUCTURE_PREDICITON  } from "./bin/subworkflows/structure.nf"
include { SEARCH                } from "./bin/subworkflows/search.nf"
include { VALIDATE              } from "./bin/subworkflows/validate.nf"
include { REPORT                } from "./bin/subworkflows/output.nf"


workflow {

    ///
    /// IMPORT & FILTERING
    ///

    // Import protein files and limit length to 1500 residues
    Channel.fromPath( params.proteins, checkIfExists: true )
    // Channel.fromPath( "./phage_data/MZ501063.1_copy/*.faa", checkIfExists: true )
    | splitFasta( record: [id:true, desc:true, seqString:true] )
    | filter { record -> record.seqString.length() < 100 }
    | set { ch_allProteins }

    FILTER( ch_allProteins )


    ///
    /// CLUSTERING
    ///

    CLUSTER( FILTER.out )


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